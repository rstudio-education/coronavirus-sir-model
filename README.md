# coronavirus-sir-model
Carl's personal by-country SIR model of the Coronavirus spread.

This model loads the latest Johns Hopkins actual data and then runs an SIR model on it to predict when the maximum infection point happens for all the countries reported by Johns Hopkins. I welcome suggestions and edits.

Directory contents are as follows:

`SIR-model.Rmd`: a parameterized static RMarkdown document that displays the input data for the world and the US and shows you the model output for the US. You can change the parameters if you want a different country.

`Shiny-SIR-Infection-model/app.R`: A Shiny app that allows you to dynamically play with all the SIR model parameters and see how they affect the forecast. You can select any country in the data set or select "All" for aggregated data for the world.


