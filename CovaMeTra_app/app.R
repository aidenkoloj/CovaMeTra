
library(tidyverse)
library(dplyr)
library(magrittr)
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

ui <- fluidPage(
    
    headerPanel("CovaMeTra"),
    
    #fileInput("file", "Upload a csv file", accept = c(".csv")),
    
    textInput(inputId = "feature",
              label = "Enter a feature name:"),
    
    tableOutput("head"),
    
    downloadButton("downloadData", "Download Table")
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



