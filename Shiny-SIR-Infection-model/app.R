#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(shiny)
library(dplyr)
library(DT)
library(ggplot2)
library(lubridate)
library(tidyr)
library(readr)
library(scales)
library(plotly)
library(glue)
library(readr)
library(RColorBrewer)
library(stringr)
library(rvest)

# Function to read the raw CSV files. The files are aggregated to the country
# level and then converted to long format

clean_jhd_to_long <- function(df) {
    df_str <- deparse(substitute(df))
    var_str <- substr(df_str, 1, str_length(df_str) - 4)
    
    df %>% group_by(`Country/Region`) %>%
        filter(`Country/Region` != "Cruise Ship") %>%
        select(-`Province/State`, -Lat, -Long) %>%
        mutate_at(vars(-group_cols()), sum) %>% 
        distinct() %>%
        ungroup() %>%
        rename(country = `Country/Region`) %>%
        pivot_longer(
            -country, 
            names_to = "date_str", 
            values_to = var_str
        ) %>%
        mutate(date = mdy(date_str)) %>%
        select(country, date, !! sym(var_str)) 
}

confirmed_raw <- read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")
deaths_raw <- read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv")

raw_coronavirus <- clean_jhd_to_long(confirmed_raw) %>%
    full_join(clean_jhd_to_long(deaths_raw)) %>%
    arrange(country, date)

coronavirus <- raw_coronavirus %>% 
    pivot_longer(confirmed:deaths, values_to = "cases", names_to = "type")

options(shiny.error = browser)

un_country_list <- read_csv("data/UNdata_Export_20200316.csv")
names(un_country_list) <- c("country", "year", "variant", "population")
un_country_list %>% filter(year == 2019)

coronavirus_countries <- coronavirus %>% 
    select(country) %>% 
    distinct(country)  ## for selecting a country

excluded_countries <- coronavirus_countries %>% 
    anti_join(un_country_list, by = c("country" = "country"))

countries_needing_translation <-read_csv("data/excluded_countries_dictionary.csv")

select_countries <- c("All", coronavirus_countries, recursive = TRUE, use.names = FALSE)
select_countries <- tibble(select_countries) %>% 
    filter(select_countries != "Taiwan*" & select_countries != "Cruise Ship")
names(select_countries) <- c("country")
select_countries <- select_countries %>% 
    arrange(country)
full_countries <- select_countries %>% 
    left_join(countries_needing_translation, by = c("country")) %>% 
    mutate(UN_country = ifelse(is.na(UN_country), country, UN_country)) %>% 
    left_join(un_country_list, by = c("UN_country" = "country")) %>% 
    filter(year == "2019") 

nytimes <- read_csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv")
nytimes_data <- nytimes %>% 
    select(-fips)
nytimes_us_cases <- nytimes_data %>% 
    pivot_longer(cases:deaths,  values_to = "cases", names_to = "type")

state_populations <- read_csv(here::here("data/us-state-populations.csv"))
select_states <- c("All", state_populations$state, recursive = TRUE, use.names = FALSE)

nytimes_us_cases <- nytimes_us_cases %>% 
    left_join(state_populations, by = "state") %>% 
    filter(!is.na(population))

run_model <- function(f, population, duration, rnaught, drate, proj_duration) {
    
    ww_model <- f
    ww_model <- ww_model %>% 
        mutate(type = ifelse(type == "death", "deaths", type))
    
    ## extend the date by project_duration days
    
    ww_wide <- ww_model %>% 
        pivot_wider(names_from= type, values_from=c(cases)) %>% 
        arrange(date)
    
    ww_date_range <- tibble(date =  ymd(seq(min(ww_wide$date), (max(ww_wide$date)), "days")))
    
    historic_range <- ww_date_range %>% 
        left_join(ww_wide, by = "date")
    
    historic_rows <- nrow(historic_range)
    
    ww_projection_dates <- tibble(date =  ymd(seq(min(ww_wide$date), (max(ww_wide$date) + proj_duration), "days")))
    
    ww_range <- ww_projection_dates %>% 
        left_join(ww_wide, by = "date")
    
    ## OK, here's the model
    ## Every day, we get new infections that's the (number of infectious people - deaths) * r0 / recovery_rate.
    ## New deaths is just the death_rate times infected people.
    ## New recoveries is the number of infected people from 14 days ago - number of deaths
    ## New number of infected people is prior number of infected people + new infections - deaths - new recoveries
    ## first, let's initialize our projected columns with the actual data.
    fnew <- ww_range %>% 
        mutate(total_infections     = infected,
               currently_infectious = infected - deaths,
               new_infections       = 0,
               new_deaths           = 0,
               total_deaths         = deaths,
               new_recoveries       = 0,
               total_recoveries     = 0,
               susceptible          = population - infected - deaths,
               ds_dt                = 0,
               di_dt                = 0,
               dr_dt                = 0) %>% 
        replace_na(list(cases_infected       = 0,
                        cases_deaths         = 0,
                        cases_recovered      = 0,
                        cumsum_infected      = 0,
                        cumsum_deaths        = 0,
                        cumsum_recovered     = 0,
                        total_infections     = 0, 
                        currently_infectious = 0,
                        new_infections       = 0, 
                        new_deaths           = 0, 
                        total_deaths         = 0,
                        new_recoveries       = 0,
                        total_recoveries     = 0)) %>% 
        arrange(date)
    
    ## SIR model parameter calculation
    
    beta <- as.numeric(rnaught / duration)
    gamma <- as.numeric(1 / duration)
    
    daily_infection_probability <- rnaught / duration
    
    for (i in (historic_rows + 1):nrow(fnew)) {
        
        ## Let's calculate the SIR model derivatives
        
        fnew$ds_dt[i] <- round(-beta * fnew$currently_infectious[i - 1] * fnew$susceptible[i - 1] / population) - fnew$new_deaths[i - 1]
        fnew$di_dt[i] <- round(-fnew$ds_dt[i] - gamma * fnew$currently_infectious[i - 1])
        fnew$dr_dt[i] <- round(gamma * fnew$currently_infectious[i - 1])
        
        ## Now update the state with derivatives
        
        fnew$new_infections[i]       <-  ifelse (fnew$di_dt[i] < 0, 0, fnew$di_dt[i])
        fnew$currently_infectious[i] <-  fnew$currently_infectious[i - 1] + fnew$di_dt[i]
        fnew$new_deaths[i]           <-  round(fnew$currently_infectious[i - 1] * drate / duration)
        fnew$new_recoveries[i]       <-  fnew$dr_dt[i]
        fnew$susceptible[i]          <-  fnew$susceptible[i - 1]      + fnew$ds_dt[i]
        fnew$total_recoveries[i]     <-  fnew$total_recoveries[i - 1] + fnew$new_recoveries[i]
        fnew$total_deaths[i]         <-  fnew$total_deaths[i - 1]     + fnew$new_deaths[i]
        fnew$total_infections[i]     <-  fnew$total_infections[i - 1] + fnew$new_infections[i]
    }
    return(fnew)
}

# Define UI for application that draws a histogram
ui <- fluidPage(

    
    # Application title
    titlePanel("Coronavirus SIM Model"),
    hr(),
    p(div(HTML("This model pulls daily data published by Johns Hopkins regarding the coronavirus infections around the world and applies a simple <em>Susceptible / Infectious / Recovery</em> or <em>SIR</em> model to that data."))), 
    p(div(HTML("Disclaimer: This simulation is for research and educational purposes only and is not intended to be a tool for decision-making."))),
    p(div(HTML('Source code <a href="https://github.com/rstudio-education/coronavirus-sir-model">is available here.</a>' ))),
    hr(),    
    
    ## Show a sidebar
    sidebarLayout(
        sidebarPanel(width = 3,
            sliderInput(inputId = "percent_susceptible_population", 
                        label = "Population Susceptible: ", 
                        min = 5, max = 100, 
                        value = 80, step = 5, post = "%"),
            sliderInput(inputId = "infection_duration", 
                        label = "Duration of Infection (Days): ", 
                        min = 2, max = 60, 
                        value = 14),
            sliderInput(inputId = "r0",
                        label = "R(0) transmissivity: ",
                        min = 0.1, max = 10,
                        value = 2.2, step = 0.1),
            sliderInput(inputId = "death_rate",
                        label = "Death Rate: ",
                        min = 0.1, max = 20, post = "%",
                        value = 1, step = 0.1),
            sliderInput(inputId = "projection_duration",
                        label = "Length of Projection (Days): ",
                        min = 10, max = 365,
                        value = 180)
        ),
        mainPanel(width = 9,
            ## Now show two tabs for outputs, one for countries and the other for US states
            tabsetPanel(
                tabPanel("Country", fluid = TRUE,
                         mainPanel(
                             selectInput(inputId = "country",
                                         label = "Select Country",
                                         choices = full_countries$country,
                                         selected = "US"),
                             h3(textOutput("country_title")),
                             h4(HTML("Next 30 Days")),
                             textOutput("caption"),
                             br(),
                             plotlyOutput("virusPlot30Days"),
                             br(),
                             h4(HTML("Peak Infectious Projection")),
                             textOutput("time_to_max"),
                             br(),
                             plotlyOutput("virusPlot"),
                         )
                ),
                tabPanel("US State", fluid = TRUE,
                         mainPanel(
                             selectInput(inputId = "state",
                                         label = "Select US State",
                                         choices = select_states,
                                         selected = "New York"),
                             h3(textOutput("state_title")),
                             h4(HTML("Next 30 Days")),
                             textOutput("state_caption"),
                             br(),
                             plotlyOutput("stateVirusPlot30Days"),
                             br(),
                             h4(HTML("Peak Infectious Projection")),
                             textOutput("state_time_to_max"),
                             br(),
                             plotlyOutput("stateVirusPlot")
                         )
                )
            )
        )
    )
)

testinputs <- function() {
    input <- tibble(country = "US",
    state = "Massachusetts",
    percent_susceptible_population = 80,
    infection_duration = 14,
    r0 = 2.2,
    death_rate = 2,
    projection_duration = 180)
    input
}

# Define server logic required to draw a histogram
server <- function(input, output) {
    derive_country_stats <- reactive({
#  derive_country_stats <- function() {     ## for debugging
        if (input$country == "All") {
            df <- coronavirus %>%
                group_by(date, type) %>% 
                summarize(cases = sum(cases, na.rm = TRUE)) %>% 
                select(date,type,cases) %>% 
                arrange(date)
        } else {
            df <- coronavirus %>%
                filter(country == input$country) %>% 
                group_by(date, type) %>%
                drop_na(cases) %>% 
                summarize(cases = sum(cases, na.rm = TRUE)) %>% 
                select(date,type,cases) %>% 
                arrange(date)
        }
        df <- df %>% 
            mutate(type = ifelse(type == "confirmed", "infected", type))
        df
    })

model_builder <- reactive({
        #   model_builder <- function(){
        
        country_info <- full_countries %>% 
            filter(country == input$country) %>% 
            head(1)
        
        total_susceptible_population <- country_info$population * input$percent_susceptible_population/100 * 1000
        df <- derive_country_stats()
        run_model(df, total_susceptible_population, input$infection_duration,
                  input$r0, input$death_rate/100, input$projection_duration)
    })
    
derive_state_stats <- reactive({     ## for debugging
## derive_state_stats <- function() {
        if (input$state == "All") {
            df <- nytimes_us_cases %>%
                group_by(date, type) %>% 
                summarize(cases = sum(cases, na.rm = TRUE)) %>% 
                select(date,type,cases) %>% 
                arrange(date)
        } else {
            df <- nytimes_us_cases %>%
                filter(state == input$state) %>% 
                group_by(date, type) %>%
                drop_na(cases) %>% 
                summarize(cases = sum(cases, na.rm = TRUE)) %>% 
                select(date,type,cases) %>% 
                arrange(date)
        }
        df <- df %>% 
            mutate(type = ifelse(type == "cases", "infected", type))
        df
}
)

 state_model_builder <- reactive({
#   state_model_builder <- function(){
        
     if (input$state == "All") {
         state_stats <- tibble(population = sum(state_populations$population))
      } else {
         state_stats <- state_populations %>% 
             filter(state == input$state) %>% 
             select(population)
     }
     state_population_susceptible <- state_stats$population * input$percent_susceptible_population/100
     df <- derive_state_stats()
     run_model(df, state_population_susceptible, input$infection_duration,
               input$r0, input$death_rate/100, input$projection_duration)
     
 }
 )

    output$virusPlot <- renderPlotly({
        theme_set(theme_minimal()) %+replace%
            theme(text=element_text(size=10,  family="Sans"))
        model_result <- model_builder()
        
        todays_values <- derive_country_stats() %>% 
            ungroup() %>% 
            filter(date == max(date)) %>% 
             mutate(type = ifelse(type == "death", "Deaths", str_to_title(type)),
                   casestring = prettyNum(abs(cases), big.mark = ","),
                   datestring = format(date, "%B %d"))
        
       todays_values_x_value <- todays_values$date - days(3)
       todays_values_string <-  paste0("Reported on ", first(todays_values$datestring), "\n", 
                                       paste(with(todays_values, glue("{type}: {casestring}\n")), collapse="\n"))
       
        final_values <- model_result %>% 
            filter(date == max(model_result$date)) %>% 
            head(n = 1)       ## in case we get multiple results
        
        ## name mapping
        
        name_map <- tribble( ~column_name, ~Category,
                             "total_infections",      "Total Infections",
                             "currently_infectious",  "Currently Infectious",
                             "new_infections",        "New Infections",
                             "new_deaths",            "New Deaths",
                             "total_deaths",          "Total Deaths",
                             "new_recoveries",        "New Recoveries",
                             "total_recoveries",      "Total Recoveries",
                             "susceptible",           "Suspectible Population")
        
        projection_long <- model_result %>% 
            pivot_longer(total_infections:susceptible)
        
        max_infectious <- model_result %>% 
            filter(model_result$currently_infectious == max(model_result$currently_infectious)) %>% 
            head(n = 1)       ## in case we get multiple results
        
        final_values <- model_result %>% 
            filter(date == max(model_result$date)) %>% 
            head(n = 1)       ## in case we get multiple results
        
        max_plot_y_value <- max(projection_long$value) * 0.90
        max_infection_text_x_value <- max_infectious$date - days(3)
        
        plot_df <- projection_long %>% 
            left_join(name_map, by = c("name" = "column_name"))
        
        g <-ggplot(plot_df, aes(x = date, y = value, color = Category, group = Category)) +
            geom_line() +
            scale_y_continuous(labels = unit_format(accuracy = 0.1, unit = "M", scale = 1e-6)) +
            scale_x_date(date_breaks = "2 weeks", limits = c(min(projection_long$date), max(projection_long$date) + days(10)), date_labels = "%B %d",
                         expand = expansion(0, 0.4)) +
            scale_color_brewer(palette = "Dark2") +
            
            ## annotate for latest reported values
            
            annotate("segment", x = todays_values$date, xend = todays_values$date, y = 0, yend = max_plot_y_value * 0.5, color = "gray20",  linetype = "dotted") +
            annotate("text",  x = todays_values_x_value, y = max_plot_y_value * 0.65, 
                     color = "gray30", size = 2.75,  hjust = "right", vjust = "inward",
                     label = todays_values_string) +

            ## and for peak infectious population
            
            annotate("segment", x = max_infectious$date, xend = max_infectious$date, y = 0, yend = max_plot_y_value * 0.9, color = "gray20",  linetype = "dashed") +
            annotate("text",  x = max_infection_text_x_value, y = max_plot_y_value, 
                     color = "gray20", size = 2.75,  hjust = "right", vjust = "inward",
                     label = paste0("Peak infectious population\n",
                                    prettyNum(round(max_infectious$currently_infectious), big.mark = ","),
                                    "\non ", format(max_infectious$date, "%B %d"))) +

            ## and for total deaths            

            annotate("text", x = max_infection_text_x_value + days(30), y = final_values$total_deaths * 3, 
                     color = "gray20", size = 2.75, hjust = "inward", vjust = "inward",
                     label = paste0("Total deaths\n", prettyNum(round(final_values$total_deaths), big.mark = ","))) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom") +
            labs(x = "", y = "", color = "")
        fig <- ggplotly(g, tooltip = c("x", "y", "group"))
        t <- list(
            family = "arial",
            size = 12)
        fig <- fig %>% layout(font=t)
    })
    
    output$virusPlot30Days <- renderPlotly({
        theme_set(theme_minimal()) %+replace%
            theme(text=element_text(size=12,  family="Sans"))
        todays_values <- derive_country_stats() %>% 
            ungroup() %>% 
            filter(date == max(date)) %>% 
            mutate(type = ifelse(type == "death", "Deaths", str_to_title(type)),
                   casestring = prettyNum(abs(cases), big.mark = ","),
                   datestring = format(date, "%B %d"))
        todays_values_x_value <- todays_values$date - days(3)
        todays_values_string <-  paste0("Reported on ", first(todays_values$datestring), "\n", 
                                        paste(with(todays_values, glue("{type}: {casestring}\n")), collapse="\n"))

        model_result <- model_builder() %>% 
            filter(date <= head(todays_values_x_value,1) + days(33))
        
        final_values <- model_result %>% 
            filter(date == max(model_result$date)) %>% 
            head(n = 1)       ## in case we get multiple results
        
        ## name mapping
        
        name_map <- tribble( ~column_name, ~Category,
                             "total_infections",      "Total Infections",
                             "currently_infectious",  "Currently Infectious",
                             "new_infections",        "New Infections",
                             "new_deaths",            "New Deaths",
                             "total_deaths",          "Total Deaths",
                             "new_recoveries",        "New Recoveries",
                             "total_recoveries",      "Total Recoveries",
                             "susceptible",           "Suspectible Population")
        
        projection_long <- model_result %>% 
            pivot_longer(total_infections:total_recoveries)
        
        max_infectious <- model_result %>% 
            filter(model_result$currently_infectious == max(model_result$currently_infectious)) %>% 
            head(n = 1)       ## in case we get multiple results
        
        final_values <- model_result %>% 
            filter(date == max(model_result$date)) %>% 
            head(n = 1)       ## in case we get multiple results
        
        max_plot_y_value <- max(projection_long$value) * 0.90
        max_infection_text_x_value <- max_infectious$date - days(3)
        
        plot_df <- projection_long %>% 
            left_join(name_map, by = c("name" = "column_name")) %>% 
            select(date, value, Category)
        
        g <-ggplot(plot_df, aes(x = date, y = value, color = Category, group = Category)) +
            geom_line() +
            scale_y_continuous(label = comma) +
            scale_x_date(date_breaks = "2 weeks", limits = c(min(projection_long$date), todays_values_x_value + days(33)), date_labels = "%B %d",
                         expand = expansion(0, 0.4)) +
            scale_color_brewer(palette = "Dark2") +
            
            ## annotate for latest reported values
            
            annotate("segment", x = todays_values$date, xend = todays_values$date, y = 0, yend = max_plot_y_value * 0.5, color = "gray20",  linetype = "dotted") +
            annotate("text",  x = todays_values_x_value, y = max_plot_y_value * 0.65, 
                     color = "gray30", size = 2.75,  hjust = "right", vjust = "inward",
                     label = todays_values_string) +
            
            ## and for peak infectious population
            
            annotate("segment", x = max_infectious$date, xend = max_infectious$date, y = 0, yend = max_plot_y_value * 0.9, color = "gray20",  linetype = "dashed") +
            annotate("text",  x = max_infection_text_x_value - days(14), y = max_plot_y_value, 
                     color = "gray20", size = 2.75,  hjust = "inward", vjust = "inward",
                     label = paste0("On ",format(max_infectious$date, "%B %d"),
                                    "\nInfected: ", prettyNum(round(max_infectious$currently_infectious), 
                                                              big.mark = ","),
                                    "\nDeaths: ", prettyNum(round(max_infectious$total_deaths), big.mark = ","))) +
            
            ## and for total deaths            
            
            annotate("text", x = max_infection_text_x_value + days(30), y = final_values$total_deaths * 3, 
                     color = "gray20", size = 2.75, hjust = "inward", vjust = "inward",
                     label = paste0("Total deaths\n", prettyNum(round(final_values$total_deaths), big.mark = ","))) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom") +
            labs(x = "", y = "", color = "")
        fig <- ggplotly(g, tooltip = c("x", "y", "group"))
        t <- list(
            family = "arial",
            size = 12)
        fig <- fig %>% layout(font=t)
        })
   
    output$stateVirusPlot <- renderPlotly({
        theme_set(theme_minimal()) %+replace%
            theme(text=element_text(size=10,  family="Sans"))
        model_result <- state_model_builder()
        
        todays_values <- derive_state_stats() %>% 
            ungroup() %>% 
            filter(date == max(date)) %>% 
            mutate(type = ifelse(type == "death", "Deaths", str_to_title(type)),
                   casestring = prettyNum(abs(cases), big.mark = ","),
                   datestring = format(date, "%B %d"))
        
        todays_values_x_value <- todays_values$date - days(3)
        todays_values_string <-  paste0("Reported on ", first(todays_values$datestring), "\n", 
                                        paste(with(todays_values, glue("{type}: {casestring}\n")), collapse="\n"))
        
        final_values <- model_result %>% 
            filter(date == max(model_result$date)) %>% 
            head(n = 1)       ## in case we get multiple results
        
        ## name mapping
        
        name_map <- tribble( ~column_name, ~Category,
                             "total_infections",      "Total Infections",
                             "currently_infectious",  "Currently Infectious",
                             "new_infections",        "New Infections",
                             "new_deaths",            "New Deaths",
                             "total_deaths",          "Total Deaths",
                             "new_recoveries",        "New Recoveries",
                             "total_recoveries",      "Total Recoveries",
                             "susceptible",           "Suspectible Population")
        
        projection_long <- model_result %>% 
            pivot_longer(total_infections:susceptible)
        
        max_infectious <- model_result %>% 
            filter(model_result$currently_infectious == max(model_result$currently_infectious)) %>% 
            head(n = 1)       ## in case we get multiple results
        
        final_values <- model_result %>% 
            filter(date == max(model_result$date)) %>% 
            head(n = 1)       ## in case we get multiple results
        
        max_plot_y_value <- max(projection_long$value) * 0.90
        max_infection_text_x_value <- max_infectious$date - days(3)
        
        plot_df <- projection_long %>% 
            left_join(name_map, by = c("name" = "column_name"))
        
        g <-ggplot(plot_df, aes(x = date, y = value, color = Category, group = Category)) +
            geom_line() +
            scale_y_continuous(labels = unit_format(accuracy = 0.1, unit = "M", scale = 1e-6)) +
            scale_x_date(date_breaks = "2 weeks", limits = c(min(projection_long$date), max(projection_long$date) + days(10)), date_labels = "%B %d",
                         expand = expansion(0, 0.4)) +
            scale_color_brewer(palette = "Dark2") +
            
            ## annotate for latest reported values
            
            annotate("segment", x = todays_values$date, xend = todays_values$date, y = 0, yend = max_plot_y_value * 0.5, color = "gray20",  linetype = "dotted") +
            annotate("text",  x = todays_values_x_value, y = max_plot_y_value * 0.65, 
                     color = "gray30", size = 2.75,  hjust = "right", vjust = "inward",
                     label = todays_values_string) +
            
            ## and for peak infectious population
            
            annotate("segment", x = max_infectious$date, xend = max_infectious$date, y = 0, yend = max_plot_y_value * 0.9, color = "gray20",  linetype = "dashed") +
            annotate("text",  x = max_infection_text_x_value, y = max_plot_y_value, 
                     color = "gray20", size = 2.75,  hjust = "right", vjust = "inward",
                     label = paste0("Peak infectious population\n",
                                    prettyNum(round(max_infectious$currently_infectious), big.mark = ","),
                                    "\non ", format(max_infectious$date, "%B %d"))) +
            
            ## and for total deaths            
            
            annotate("text", x = max_infection_text_x_value + days(30), y = final_values$total_deaths * 3, 
                     color = "gray20", size = 2.75, hjust = "inward", vjust = "inward",
                     label = paste0("Total deaths\n", prettyNum(round(final_values$total_deaths), big.mark = ","))) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom") +
            labs(title = "State projection", x = "", y = "", color = "")
        fig <- ggplotly(g, tooltip = c("x", "y", "group"))
        t <- list(
            family = "arial",
            size = 12)
        fig <- fig %>% layout(font=t)
    })
    
    output$stateVirusPlot30Days <- renderPlotly({
        theme_set(theme_minimal()) %+replace%
            theme(text=element_text(size=10,  family="Sans"))
        
        todays_values <- derive_state_stats() %>% 
            ungroup() %>% 
            filter(date == max(date)) %>% 
            mutate(type = ifelse(type == "death", "Deaths", str_to_title(type)),
                   casestring = prettyNum(abs(cases), big.mark = ","),
                   datestring = format(date, "%B %d"))
        
        todays_values_x_value <- todays_values$date - days(3)
        todays_values_string <-  paste0("Reported on ", first(todays_values$datestring), "\n", 
                                        paste(with(todays_values, glue("{type}: {casestring}\n")), collapse="\n"))
        model_result <- state_model_builder() %>% 
            filter(date <= head(todays_values_x_value,1) + days(33))
        
        final_values <- model_result %>% 
            filter(date == max(model_result$date)) %>% 
            head(n = 1)       ## in case we get multiple results
        
        ## name mapping
        
        name_map <- tribble( ~column_name, ~Category,
                             "total_infections",      "Total Infections",
                             "currently_infectious",  "Currently Infectious",
                             "new_infections",        "New Infections",
                             "new_deaths",            "New Deaths",
                             "total_deaths",          "Total Deaths",
                             "new_recoveries",        "New Recoveries",
                             "total_recoveries",      "Total Recoveries",
                             "susceptible",           "Suspectible Population")
        
        projection_long <- model_result %>% 
            pivot_longer(total_infections:susceptible)
        
        max_infectious <- model_result %>% 
            filter(model_result$currently_infectious == max(model_result$currently_infectious)) %>% 
            head(n = 1)       ## in case we get multiple results
        
        final_values <- model_result %>% 
            filter(date == max(model_result$date)) %>% 
            head(n = 1)       ## in case we get multiple results
        
        max_plot_y_value <- max(projection_long$value) * 0.90
        max_infection_text_x_value <- max_infectious$date - days(3)
        
        plot_df <- projection_long %>% 
            left_join(name_map, by = c("name" = "column_name"))
        
        g <-ggplot(plot_df, aes(x = date, y = value, color = Category, group = Category)) +
            geom_line() +
            scale_y_continuous(labels = unit_format(accuracy = 0.1, unit = "M", scale = 1e-6)) +
            scale_x_date(date_breaks = "2 weeks", limits = c(min(projection_long$date), max(projection_long$date) + days(10)), date_labels = "%B %d",
                         expand = expansion(0, 0.4)) +
            scale_color_brewer(palette = "Dark2") +
            
            ## annotate for latest reported values
            
            annotate("segment", x = todays_values$date, xend = todays_values$date, y = 0, yend = max_plot_y_value * 0.5, color = "gray20",  linetype = "dotted") +
            annotate("text",  x = todays_values_x_value, y = max_plot_y_value * 0.65, 
                     color = "gray30", size = 2.75,  hjust = "right", vjust = "inward",
                     label = todays_values_string) +
            
            ## and for peak infectious population
            
            annotate("segment", x = max_infectious$date, xend = max_infectious$date, y = 0, yend = max_plot_y_value * 0.9, color = "gray20",  linetype = "dashed") +
            annotate("text",  x = max_infection_text_x_value, y = max_plot_y_value, 
                     color = "gray20", size = 2.75,  hjust = "right", vjust = "inward",
                     label = paste0("Peak infectious population\n",
                                    prettyNum(round(max_infectious$currently_infectious), big.mark = ","),
                                    "\non ", format(max_infectious$date, "%B %d"))) +
            
            ## and for total deaths            
            
            annotate("text", x = max_infection_text_x_value + days(30), y = final_values$total_deaths * 3, 
                     color = "gray20", size = 2.75, hjust = "inward", vjust = "inward",
                     label = paste0("Total deaths\n", prettyNum(round(final_values$total_deaths), big.mark = ","))) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom") +
            labs(title = "State projection", x = "", y = "", color = "")
        fig <- ggplotly(g, tooltip = c("x", "y", "group"))
        t <- list(
            family = "arial",
            size = 12)
        fig <- fig %>% layout(font=t)
    })

       
    output$country_title <- renderText({
        glue("{input$country} Coronavirus Forecast")
            })
    
    output$state_title <- renderText({
        glue("{input$state} Coronavirus Forecast")
    })
    
    output$time_to_max <- renderText({
        model_result <- model_builder()
        max_infectious <- model_result %>% 
            filter(model_result$currently_infectious == max(model_result$currently_infectious)) %>% 
            head(n = 1)       ## in case we get multiple results
        
        first_string <- glue("{max_infectious$date - today()} days until infections peak in {input$country}.")
        country_info <- full_countries %>% 
            filter(country == input$country) %>% 
            head(1)
        total_susceptible_population <- country_info$population * input$percent_susceptible_population/100 * 1000
        second_string <- glue("Susceptible population = {prettyNum(round(total_susceptible_population/1e6,1), big.mark = ',')}M, R0 = {input$r0}, infection duration = {input$infection_duration}, death rate = {input$death_rate}%")
        return(c(first_string, second_string))
        })

    output$state_time_to_max <- renderText({
        model_result <- state_model_builder()
        max_infectious <- model_result %>% 
            filter(model_result$currently_infectious == max(model_result$currently_infectious)) %>% 
            head(n = 1)       ## in case we get multiple results
        
        first_string <- glue("{max_infectious$date - today()} days until infections peak in {input$state}.")
        if (input$state == "All") {
            state_stats <- tibble(population = sum(state_populations$population))
        } else {
            state_stats <- state_populations %>% 
                filter(state == input$state) %>% 
                select(population)
        }
        state_population_susceptible <- state_stats$population * input$percent_susceptible_population/100
        second_string <- glue("Susceptible population = {prettyNum(round(state_population_susceptible/1e6,1), big.mark = ',')}M, R0 = {input$r0}, infection duration = {input$infection_duration}, death rate = {input$death_rate}%")
        return(c(first_string, second_string))
    })
        
    output$caption <- renderText({
        country_info <- full_countries %>% 
            filter(country == input$country) %>% 
            head(1)
        total_susceptible_population <- country_info$population * input$percent_susceptible_population/100 * 1000
        glue("Susceptible population = {prettyNum(round(total_susceptible_population/1e6,1), big.mark = ',')}M, R0 = {input$r0}, infection duration = {input$infection_duration}, death rate = {input$death_rate}%")
    })

    output$state_caption <- renderText({
        if (input$state == "All") {
            state_stats <- tibble(population = sum(state_populations$population))
            state_population_susceptible <- state_stats$population * input$percent_susceptible_population/100
            glue("US susceptible population = {prettyNum(round(state_population_susceptible/1e6,1), big.mark = ',')}M, R0 = {input$r0}, infection duration = {input$infection_duration}, death rate = {input$death_rate}%")
        } else {
            state_stats <- state_populations %>% 
                filter(state == input$state) %>% 
                select(population)
            state_population_susceptible <- state_stats$population * input$percent_susceptible_population/100
            glue("{input$state} Susceptible population = {prettyNum(round(state_population_susceptible/1e6,1), big.mark = ',')}M, R0 = {input$r0}, infection duration = {input$infection_duration}, death rate = {input$death_rate}%")
        }
    })
    
    }

# Run the application 
shinyApp(ui = ui, server = server)
