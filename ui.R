######################################################################
###Author: aekosel
###Purpose: This is the UI code for R Shiny. It also requires
###the server code for use.
######################################################################

####################
###load libraries###
####################
require(shiny)
require(KernSmooth)
dat = read.csv("simulated_data.csv")


####################
###user interface###
####################      
shinyUI(pageWithSidebar(

	#Application title
	headerPanel("Localized Patient Outcome Predictions"),

	sidebarPanel(
 
		#choose outcome, x1, x2 variables
		selectInput("myout", "Choose Y Variable", names(dat)),
		selectInput("myvar1", "Choose X1 Variable", names(dat)),
		selectInput("myvar2", "Choose X2 Variable", names(dat)),

		#scaled?
		selectInput("scaled", "Scale covariates?", c("No", "Globally", "Locally")),

		#select point of interest
		numericInput("my.point.x1", "My Point X1", 0),
		numericInput("my.point.x2", "My Point X2", 0),

		#sliders for x1, x2 plot limits
		uiOutput("x1.lim"),
		uiOutput("x2.lim"),

		#slider for points in neighborhood
		uiOutput("n"),
 
		#slider for histogram breaks
		sliderInput("histboxes", "No. of Breaks in Histogram:", 
			min=1, max=50, value=20),

		#slider for histogram quantile
		sliderInput("histquant", "Quantile in Histogram:", 
			min=0, max=100, value=50)

	),

	mainPanel(
		#display covariates
		div(class="span6", plotOutput("myplot")),

		#display outcome
		div(class="span6", plotOutput("myhist")),

		#display uncertainty
		div(class="span3", tableOutput("citable")),

		#display quantiles
		div(class="span3", tableOutput("quantiletable")),

		#display neighborhood edges
		div(class="span3", tableOutput("cornertable")),

		#display number of neighborhood points
		div(class="span3", tableOutput("inboxtable")),

		#display scaling
		div(class="span4", tableOutput("boxsizetable"))

	)
))

