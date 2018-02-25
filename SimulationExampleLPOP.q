#
# Example of steps/commands needed to execute LPOP shiny app
#
###### [1] -- set the desired data to "lpop.data"
#
lpop.data <- read.csv("simulated_data.csv")
#
##### [2] -- load the server and UI
#
source("server.R")
source("ui.R")
#
##### [3] -- launch the app!
#
runApp()
#
#end-of-file...
