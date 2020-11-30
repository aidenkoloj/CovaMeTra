
library(tidyverse)
library(dplyr)
library(reshape2)
library(magrittr)
library(shinycssloaders)
library(shinythemes)
library(shiny)
library(plotly)
options(shiny.maxRequestSize = 600*1024^2)

rna_seq = read_csv("NIDDK_RNA_seq_with_replicate_ID.csv")
metabolites = melt(read.csv("NIDDK_soup_neg_11_02.csv"))
metabolites = metabolites %>% rename(feature = variable, abundance = value, replicate = Replicate) %>% select(replicate, feature, abundance)
metabolites$feature <- paste(metabolites$feature, sep="_","metabolite")
rna_seq = rna_seq %>% rename(feature = ext_gene, abundance = normalized_counts) %>% select(replicate, feature, abundance)

metra_data <- rbind(metabolites, rna_seq)


covametra = function(x){
    
    input_data <- 
        metra_data %>% 
        #Select columns of interest
        filter(feature == x) 
    
    df <- data.frame(a = c(input_data$feature), b = c(input_data$abundance)) 
    df <- bind_rows(replicate((nrow(metra_data)/24), df, simplify = FALSE))
    combined <- cbind(metra_data, df) 
    combined_r2 <- combined %>% rename(input_names = a, input_abundance = b) %>% group_by(feature) %>% 
        summarise(r_2 = cor(abundance,input_abundance), method = "pearson") 
    combined_r2$r_2 = combined_r2$r_2*combined_r2$r_2
    combined_r2 <- combined_r2 %>% 
        arrange(desc(r_2)) 
    return(combined_r2)
    
    
}


covametra_ratio = function(x, y){
    
    numerator <- 
        metra_data %>% 
        #Select columns of interest
        filter(feature == x) %>%
        rename(featurex = feature, abundancex = abundance) %>%
        select(featurex,abundancex)
    
    denominator <-  
        metra_data %>% 
        #Select columns of interest
        filter(feature == y) %>%
        rename(featurey = feature, abundancey = abundance) %>%
        select(featurey,abundancey)
    
    input_data <- cbind(numerator, denominator) %>% mutate(ratio = (abundancex/abundancey))
    input_data$f1_f2 = paste(input_data$featurex,input_data$featurey)
    #print(input_data)
    
    df <- data.frame(a = c(input_data$f1_f2), b = c(input_data$ratio)) 
    df <- bind_rows(replicate((nrow(metra_data)/24), df, simplify = FALSE))
    combined <- cbind(metra_data, df) 
    combined_r2 <- combined %>% rename(input_names = a, input_abundance = b) %>% group_by(feature) %>% 
        summarise(r_2 = cor(abundance,input_abundance), method = "pearson") 
    combined_r2$r_2 = combined_r2$r_2*combined_r2$r_2
    combined_r2 <- combined_r2 %>% 
        arrange(desc(r_2)) 
    return(combined_r2)
    
    
}


cbPalette <- c(rep("#cf6090",6), rep("#3853a4",6), rep("#78c6cf",6), rep("#faa41a",6))

plot_2features = function(f1, f2){
    
    title = paste(f1, f2, sep ="/")
    
    df1 <- metra_data %>% filter(feature == f1) %>% rename( df1_value = abundance)
    df2 <- metra_data %>% filter(feature == f2) %>% rename( df2_value = abundance) %>% select(df2_value)
    df_combined <- cbind(df1, df2) 
    print(df_combined)
    contains_r2_value = lm(df1_value ~ df2_value, data=df_combined)
    
    ggplot(df_combined, aes(df1_value, df2_value))+
        geom_point(aes(color = replicate), size = 2)+
        #ggtitle(title)+
        theme(plot.title = element_text(hjust = 0.5, size =20, vjust = 5, face = "bold"))+
        theme_bw() + theme(panel.border = element_blank(),
                           text=element_text(size=10,  family="Arial"), legend.position = "none")+
        theme(plot.title = element_text(hjust = 0.5, size =10, vjust = 5, face = "bold"))+
        theme(plot.margin=unit(c(1,1.5,1.5,1.2),"cm"))+
        geom_smooth(method = "lm", col = "gray", se = F)+
        labs( x = as.character(f1), y = as.character(f2))+
        ### PROBLEM DISPALYING R2 is plotly geom_text(x = 0, y = 0,hjust = 0, vjust = 1, label = paste("R2 = ",signif(summary(contains_r2_value)$r.squared, 5)))+
        scale_colour_manual(values=cbPalette)
}

plot1feature <- function(x){
    df  <- metra_data %>% filter(feature == x)
    ggplot(df, aes(replicate,abundance))+
        geom_point(aes(color = replicate), size = 2)+
        ggtitle(x)+  theme(plot.title = element_text(hjust = 0.5, size =20, vjust = 5, face = "bold"))+
        theme_bw() + theme(panel.border = element_blank())+
        theme(plot.title = element_text(hjust = 0.5, size =20, vjust = 5, face = "bold"))+
        theme(plot.margin=unit(c(1,1.5,1.5,1.2),"cm"))+
        scale_colour_manual(values=cbPalette)
}

ui <- fluidPage(
    
    navbarPage("CovaMeTra", theme = shinytheme("lumen"),
               
               tabPanel("Correlate Features", fluid = TRUE, icon = icon("connectdevelop"),
                        textInput(inputId = "feature",
                                  label = "Perform a correlation, enter a feature name:"),
                        tableOutput("r2"),
                        downloadButton("downloadData1", "Download Table")),
               
               tabPanel("Correlate a Feature Ratio", fluid = TRUE, icon = icon("connectdevelop"),
           #Add text saying "create ratio"
                        textInput(inputId = "numerator",
                                  label = "Enter a feature numerator:"),
                        textInput(inputId = "denominator",
                                  label = "Enter a feature denominator"),
           ## Insert an action button
                        actionButton("calc_ratio", "Correlate Ratio"),
                        tableOutput("r2_ratio"),
                        downloadButton("downloadData2", "Download Table")),
               

                tabPanel("Plot Two Features", fluid = TRUE, icon = icon("grip-lines"),
                        textInput(inputId = "f1",
                                    label = "Enter a feature name (f1):"),
                        textInput(inputId = "f2",
                                  label = "Enter a feature name (f2):"),
                        mainPanel(plotlyOutput(outputId = "features_graph"), width = "100%")),
           
                tabPanel("Plot Single Feature", fluid = TRUE, icon = icon("grip-lines"),
                        textInput(inputId = "fsingle",
                                label = "Enter a feature"),
                        mainPanel(plotOutput(outputId = "fsingle_graph"), width = "100%")),
                        
                        
                tabPanel("Feature search", fluid = TRUE, icon = icon("search"),
                        textInput(inputId = "fsearch",
                                  label ="Search for a feature: "),
                        mainPanel(tableOutput("search_data"), width = "100%")),
           
                tabPanel("About", fluid = TRUE,
                    fluidRow(
                        column(6,
                               #br(),
                               h4(p("About CovaMeTra")),
                               h5(p("Correlation of transcriptomic data for the generation of gene regulatory networks is a well established method.
                                    On the other hand, the incorporation of metabolomic data in these networks has not been fully realized. 
                                    This program attempts to bridge these two datasets, creating a model of gene expression which is reflected in the measured metabolites.
                                    Combining these two datasets provides a tool with which to create dynamic metabolomic/transciptomic networks.")),
                               br(),
                               h4(p("About the data")),
                               h5(p("The data for this study consists of 4 strains of C. elegans: BRC20067, CB4856, DL238, and N2. 
                                    In graphs, the strains are colored Salmon, Dark Blue, Light Blue, and Yellow, repectively.
                                    Six biological replicates were grown for each strain. Worms were grown in liquid culture and fed lyophilized HB101 bacteria. 
                                    Metabolomic data was obtained for the supernatant and pelleted material from the cultures and data was collected in positive
                                    and neagtive mode. MS2 fragmentation data is also available. For more information on the dataset contact ark259@cornell.edu."),
                                  p("The source code for this Shiny app and the data is available ", a("on github", href = "https://github.com/aidenkoloj/CovaMeTra"), "."))
                               
                               #hr(),
                               
                        ),
                        column(6,
                               h4(p("About the Author")),
                               h5(p("Aiden Kolodziej is a senior studying biological sciences at Cornell University. He has spent the last year in the Schroeder Lab at the Boyce Thompson Institute learning how to obtain and anaylze metabolomic data. He is interested in combining genomic, transcriptomic, and metabolomic datasets for the purpose of developing more detailed biological system models"),
                               ),
                               HTML('<img src="DSC06163.jpg", height ="200px"'),
                               br()
                        )
                    ),
                    br(),
                    hr(),
                    h5("Sources:"),
                    h6(
                        p("Project a joint collaboration between the Andersen Lab (Northwestern), Walhout Lab (UMASS Medical), and Schroeder Lab (Cornell) "),
                        p("Shiny theme/layout from ",
                          a("gpilgrim2670", href = "https://shiny.rstudio.com/gallery/ncaa-swim-team-finder.html"),
                        )),
                        
                        

)
)
)

    
server <- function(input, output, session) {
    
    feature_name_input <- reactive({
        if(is.null(input$feature)){
            return(NULL)
        }
        input$feature
    })
    
    
    r_2 <- reactive({    
        req(input$feature)
        covametra(feature_name_input())
        
        
        
    })
    
    output$r2 <- renderTable({
        r_2()
    })
    
    observeEvent(input$calc_ratio,{    
        r_2_ratio <- 
        covametra_ratio(input$numerator, input$denominator)
        output$r2_ratio <- renderTable({
            r_2_ratio
        })
    } )


    
    output$search_data <- renderTable({
        metra_data %>% filter(feature == input$fsearch)
    })
    
    output$features_graph <- renderPlotly({ 
        width = 200
        height = 200 
        req(input$f1)
        req(input$f2)
        ggplotly(plot_2features(input$f1,input$f2))
        
    })
    
    
    output$fsingle_graph <- renderPlot({ 
        width = 50
        height = 50
        req(input$fsingle)
        plot1feature(input$fsingle)

    })
    
    output$downloadData1 <- downloadHandler(
        filename = function() {
            paste(input$feature, ".csv", sep = "")
        },
        content = function(file) {
            write.csv(r_2(), file, row.names = FALSE)
        }
    )
    
    output$downloadData2 <- downloadHandler(
        filename = function() {
            paste(input$feature_ratio, ".csv", sep = "")
        },
        content = function(file) {
            write.csv(XXX, file, row.names = FALSE)
        }
    )
}

# combine UI and server ---------------------------------------------------

shinyApp(ui,server)



