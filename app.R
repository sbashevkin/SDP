require(leaflet)
require(leaflet.minicharts)
require(readxl)
require(lubridate)
require(tidyverse)
require(shiny)
require(shinyWidgets)

ui<-fluidPage(
  titlePanel("Bivalve visualizer"),
  fluidRow(column(5, radioButtons("Interval", "Date interval", choices=c("Month", "Year"), selected="Month")),
           column(7, uiOutput("select_Date"))),
  fluidRow(column(12, leafletOutput("Mapplot", width = "100%", height = "100%"))),
  tags$head(tags$style("#Mapplot{height:80vh !important;}"))
)

server <- function(input, output, session) {
  
  Stations<-read_excel("Bivalves app/1975-18 CPUE bivalves only, 2019Sept9.xlsx",
                       sheet = "75-17 station locations", skip=1)%>%
    select(Station=Site_Code, Latitude, Longitude)
  
  Biv<-read_excel("Bivalves app/1975-18 CPUE bivalves only, 2019Sept9.xlsx",
                  sheet = "75-18 CPUE per m2", skip=1)%>%
    select(Date, Station=StationCode, `Potamocorbula amurensis`, `Corbicula fluminea`)%>%
    gather(key="Taxa", value="CPUE", -Station, -Date)%>%
    mutate(Year=floor_date(Date, unit = "year"),
           MonthYear=floor_date(Date, unit = "month"))%>%
    left_join(Stations, by="Station")
  
  Bivsum<-reactive({
    Biv%>%
      {
        if(input$Interval=="Year"){
          group_by(., Year, Taxa, Latitude, Longitude)
        } else{
          group_by(., MonthYear, Year, Taxa, Latitude, Longitude)
        }
      }%>%
      summarise(CPUE=mean(CPUE, na.rm=T))%>%
      ungroup()
  })
  
  output$select_Date <- renderUI({
    
    choice_Date <- reactive({
      Bivsum()%>%
        {
          if(input$Interval=="Year"){
            pull(., Year)
          } else{
            pull(., MonthYear)
          }
        }%>%
      unique()
    })
    
    sliderTextInput("Date",
                "Select date:",
                choices=choice_Date(),
                animate=animationOptions(interval=ifelse(input$Interval=="Year", 1000, 200)), width="100%")
  })
  
  mapdata<-reactive({
    Bivsum()%>%
      spread(key=Taxa, value=CPUE)%>%
      mutate(Total=`Corbicula fluminea`+`Potamocorbula amurensis`)
  })
  
  filteredmapdata<-reactive({
    mapdata()%>%
      {
        if(input$Interval=="Year"){
          filter(., Year==as_date(input$Date))
        } else{
          filter(., MonthYear==as_date(input$Date))
        }
      }
  })
  
  pal<-reactive({
    colorFactor(c("#1b9e77", "#7570b3"), unique(Bivsum()$Taxa))
  })
  
  Mapplot<-reactive({
    leaflet(data = mapdata())%>%
      addProviderTiles("Esri.WorldGrayCanvas")%>%
      fitBounds(~min(Longitude, na.rm=T), ~min(Latitude, na.rm=T), ~max(Longitude, na.rm=T), ~max(Latitude, na.rm=T))%>%
      addLegend("topleft", pal = pal(), values = unique(Bivsum()$Taxa))
  })
  
  output$Mapplot <- renderLeaflet({
    Mapplot()
  })
  
  observeEvent(input$Date, {
    req(nrow(filteredmapdata())>1)
    colors <- c("#1b9e77", "#7570b3")
    zeros<-filter(filteredmapdata(), Total==0)
    map<-leafletProxy("Mapplot", session, data = filteredmapdata())%>%
      clearMinicharts() %>%
      clearMarkers()%>%
      addMinicharts(lng = filteredmapdata()$Longitude, lat = filteredmapdata()$Latitude,
                    type = "pie",
                    chartdata = filteredmapdata()%>%select_at(vars(unique(Bivsum()$Taxa)))%>%as.matrix(), 
                    colorPalette = colors, transitionTime = 0, opacity=0.5, width=60*(sqrt(filteredmapdata()$Total)/sqrt(max(mapdata()$Total))), legend=F)%>%
      addCircleMarkers(lng = zeros$Longitude, lat = zeros$Latitude, fillColor="Black", radius=3, stroke=0, fillOpacity = 1)
  })
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)