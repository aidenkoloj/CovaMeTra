
library(tidyverse)
library(dplyr)
library(magrittr)
library(shinycssloaders)
library(shinythemes)
library(shiny)
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
        geom_text(x = -Inf, y = Inf,hjust = 0, vjust = 1, label = paste("R2 = ",signif(summary(contains_r2_value)$r.squared, 5)))+
        scale_colour_manual(values=cbPalette)
}

ui <- fluidPage(
    
    navbarPage("CovaMeTra", theme = shinytheme("lumen"),
               
               tabPanel("Correlate Features", fluid = TRUE, icon = icon("connectdevelop"),
                        textInput(inputId = "feature",
                                  label = "Enter a feature name:"),
                                  tableOutput("head"),
                                  downloadButton("downloadData", "Download Table")),
               

                tabPanel("Plot Two Features", fluid = TRUE, icon = icon("grip-lines"),
                        textInput(inputId = "f1",
                                    label = "Enter a feature name (f1):"),
                        textInput(inputId = "f2",
                                  label = "Enter a feature name (f2):"),
                mainPanel(
                    
                    plotOutput(outputId = "features_graph"), width = "100%")
                
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
    
    output$head <- renderTable({
        r_2()
    })
    
    output$features_graph <- renderPlot({ 
        width = 200
        height = 200 
        req(input$f1)
        req(input$f2)
        plot_2features(input$f1,input$f2)
        
    })
    
    output$downloadData <- downloadHandler(
        filename = function() {
            paste(input$gene, ".csv", sep = "")
        },
        content = function(file) {
            write.csv(r_2(), file, row.names = FALSE)
        }
    )
}

# combine UI and server ---------------------------------------------------

shinyApp(ui,server)



