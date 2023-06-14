### Shiny App ###

#install.packages("shiny")
library("shiny")

## Does this need to be run every time?
#renv::restore()
#renv::activate()


#--------------#
#     Prep     #
#--------------#

source("code/02_generate_lotteries/helper.R")
source("code/01_mle_model/helpers/complexity_feature_functions.R")


#------------#
#     UI     #
#------------#

ui <- pageWithSidebar(
  # App title ----
  headerPanel("App"),
  # Sidebar panel for inputs ----
  sidebarPanel(   
    numericInput("head_rows", "Rows", value = 5, min = 1, max = 10, step = 1),
  ),
  # Main panel for displaying outputs ----
  mainPanel(
    fileInput("upload", "Upload file", 
              buttonLabel = "Upload files...", 
              accept = ".csv"),
    #numericInput("n", "maximum number of states", value = 7, min = 1, step = 1),
    tableOutput("preview"),
    downloadButton("download", "Download .csv")
    
    # TO DO: MAKE DONWLOAD CUSTOMIZABLE: where to save etc.
  )
)

#------------#
#   SERVER   #
#------------#
server <- function(input, output) {
  
  # LOAD DATA & PERFORM VALIDATIONS  -------------------------------------------------------------
  data_raw <- reactive({
    
    # Validate that a data set has been uploaded
    # Check this StackOverflow response to only display validation once: https://stackoverflow.com/questions/56041113/why-does-shiny-display-validation-message-twice
    shiny::validate(
      need(input$upload != "", "Please upload a data set.")
    )
    
    # Read in the data
    shiny::req(input$upload)
    file <- read.csv(input$upload$datapath, header = TRUE)
    if (is.null(file)) return(NULL)
    
    # Validate the data set's columns
    reqs = c(glue("p_a_{1:7}"), glue("p_b_{1:7}"), glue("x_a_{1:7}"), glue("x_b_{1:7}"))
    names = colnames(file)
    shiny::validate(
      need(all(reqs %in% names), "Please ensure the file contains the right columns.")
    )
    
    # shiny::validate(
    #   need(is.data.frame(file), "File cannot be loaded as a data frame.")
    # )
    
    file
  })
  
  # Execute analysis steps within reactive frame using functions
  # 1) CREATE COMPLEXITY FEATURES -------------------------------------------------------------
  data_complexity_features <- reactive({
    df = data_raw()
    
    #  1.1 Add problem column - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    df$Problem = 1:nrow(df)
    df = df %>% relocate("Problem", .before = 1)
    
    df_primitive <- df
    
    #  1.2 Repack p/x columns - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    types <- c("p", "x")
    lots <- c("a", "b")
    for (type in types) {
      for (lot in lots) {
        repacked <- glue ("{type}_{lot}")
        name_scheme <- glue("{repacked}_")
        df <- repack_column(repacked, name_scheme, df)
        for(i in 1:nrow(df)) {
          li_val <- df[[repacked]][i] %>% unlist
          if(all(is.na(li_val))) {
            df[[repacked]][i] <- NA
          } else {
            df[[repacked]][i] <- list(li_val[!is.na(li_val)])
          }
        }
      }
    }
    
    #  1.3 Generate correlated state variables - - - - - - - - - - - - - - - - - - 
    x_a_cor <- c()
    x_b_cor <- c()
    p_ab_cor <- c()
    nstates_cor <- c()
    
    for (i in seq(nrow(df))) {
      x_a_i <- df$x_a[i] %>% unlist
      x_b_i <- df$x_b[i] %>% unlist
      p_a_i <- df$p_a[i] %>% unlist
      p_b_i <- df$p_b[i] %>% unlist
      
      cor_df <- correlate_states(x_a_i, x_b_i, p_a_i, p_b_i)
      x_a_cor_i <- cor_df$pay_a
      x_b_cor_i <- cor_df$pay_b
      p_ab_cor_i <- cor_df$p_ab_cor
      nstates_cor_i <- nrow(cor_df)
      
      x_a_cor <- c(x_a_cor, list(x_a_cor_i))
      x_b_cor <- c(x_b_cor, list(x_b_cor_i))
      p_ab_cor <- c(p_ab_cor, list(p_ab_cor_i))
      nstates_cor <- c(nstates_cor, nstates_cor_i)
    }
    df$x_a_cor <- x_a_cor
    df$x_b_cor <- x_b_cor
    df$p_ab_cor <- p_ab_cor
    df$nstates_cor <- nstates_cor
    
    #  1.4 Construct features - - - - - - - - - - - - - - - - - - - - - - - - - - 
    old_names <- names(df)
    df <- build_features(df)
    
    feat_names <- names(df) %>% setdiff(old_names)
    df <- df %>% select(all_of(c("Problem", feat_names)))
    
    
    #  1.5 Store the data in data frame - - - - - - - - - - - - - - - - - - - - - - 
    final_df_complexity_features <- df_primitive %>% left_join(df)
    
    final_df_complexity_features
    #res2 = df_to_save %>% tidyr::unnest_wider(col = colnames(df_to_save), names_sep = "__")
    #res2 = df_to_save %>% select(-c("x_a_cor", "x_b_cor", "p_cor"))
    
  })
  
  
  # 2) CREATE COMPLEXITY INDICES -------------------------------------------------------------
  data_complexity_indices <- reactive({
    
    #  2.1 Load data - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    outcome_var <- "inaccurate_obj"
    
    # Load pre-processed input data
    df = data_complexity_features()
    
    # Load features
    x_train_ave_lasso <- readRDS(here::here("data", "lasso_features_ave_train_df.rds"))
    
    # Load model
    load_fp_ave <- save_fp <- here::here(
      "data", "models",
      glue("{outcome_var}__LASSO__ave.rds")
    )
    model_ave <- readRDS(load_fp_ave)
    
    #  2.2 Get LASSO "Complexity Index" - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    # AVE
    drop_vars <- setdiff(names(df), names(x_train_ave_lasso))
    x_lasso_ave <- df %>%
      select(-all_of(drop_vars))
    
    # manual fill some columns with 0 that are not defined in x_lasso_ave, Cassidy's code still being updated hence the inconsistencies
    x_lasso_ave$smry_both_mixed_nondom__ab = 0
    x_lasso_ave$smry_any_mixed_nondom__ab = 0
    x_lasso_ave$smry_two_types_dom__ab = 0
    
    # order columns to match  ## -- TO DO: DOUBLE CHECK 
    x_lasso_ave <- x_lasso_ave[, names(x_train_ave_lasso)] %>%
      as.matrix
    #ci_ave <- predict(model_ave, x_lasso_ave, s = model_ave$lambda.1se) ### WHEN CALLING PREDICT IT IS CALLED FROM STATS PACKAGE, NO? IS that a problem?
    #ci_ave <- glmnet::predict.glmnet(model_ave, x_lasso_ave, s = "lambda.1se") ## This throws an error
    library(glmnet)
    ci_ave <- predict(model_ave, x_lasso_ave, s = "lambda.1se") 
    
    # SAVE
    # xwalk
    ci_xwalk <- data.frame(
      Problem = df$Problem,
      drop = TRUE
    )
    ave_pred_name <- glue("ci_ave_{outcome_var}")

    ci_xwalk[[ave_pred_name]] <- ci_ave[,1]

    ci_xwalk <- ci_xwalk %>%
      select(-drop)
    
    final_df_complexity_indices = ci_xwalk
    final_df_complexity_indices
  })
  
  
  

  # Output the head of the dataset -- commented out for now, because list cells do not display and stop the App from running
  output$preview <- renderTable({
    if (!is.null(data_complexity_indices())){
      head(data_complexity_indices(), input$head_rows)
    }
  })
  
  # Download button
  output$download <- downloadHandler(
    filename = function() {
      paste0(input$upload, "_FINAL.csv")
    },
    content = function(file) {
      write.csv(data_complexity_indices(), file)
    }
  )
  
}

shinyApp(ui, server)

#renv::deactivate()
