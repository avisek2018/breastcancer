#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(ggplot2)
library(plotly)
library(tidyverse)
library(caret)
library(reshape2)
library(ggbiplot)
library(DT)
library(ModelMetrics)
library(class)

source("source.R")


server <- function(input, output, session) {
    
    # # Simulate work being done for 1 second
    # Sys.sleep(5)
    # 
    # # Hide the loading message when the rest of the server function has executed
    # hide(id = "loading-content", anim = TRUE, animType = "fade")   
    # show("app-content")
    
    newVar <- reactive({
        input$predictor
        input$diagram
        input$density
        input$position
        input$pcs1
        input$pcs2
        
        
    })
    
    
    #create plot
    output$plot1 <- renderPlotly({
        xValue <- input$predictor
        #print(xValue)
        if (input$diagram == "Bar Chart"){
            g <- ggplot(data = data, aes(x = pull(data, xValue)))
            g <- g + geom_bar(aes(fill = Class), position = input$position) + 
                labs(x = xValue, title = paste0("Bar Chart for ", xValue)) + 
                scale_fill_discrete(labels = c("Non Malignant", "Malignant"))
            ggplotly(g)
        } else if (input$diagram == "Histogram"){
            g <- ggplot(data = data, aes(x = pull(data, xValue))) + labs( x = xValue, title = paste0("Histogram for ", xValue))
            if (input$density) {
                g <- g + geom_histogram(bins = 20, aes(y = ..density.., fill = Class), position = input$position) + 
                    geom_density(adjust = 0.25, alpha = 0.5, aes(fill=Class), position = input$position)
                ggplotly(g)
            }
            else { g + geom_histogram(bins = 20, aes(y = ..density.., fill = Class), position = input$position) }
        }else if(input$diagram == "BoxPlot"){
            g <- ggplot(data = data, aes(x = Class, y= pull(data, xValue))) 
            g <- g + geom_boxplot() + geom_jitter(aes(color = Class)) + labs(y = xValue, title = paste0("Boxplot for ", xValue))
            ggplotly(g)
        }
    })
    
    output$summary <- renderPrint({
        summary(data[input$predictor])
    })
    
    # output$allboxplot <- renderPlot({
    #     df.m <- melt(data, id.vars = "Class")
    #     df.m <- df.m %>% drop_na()
    #     df.m$Class <- as.factor(df.m$Class)
    #     g <- ggplot(data = df.m, aes(x=variable, y=value)) 
    #     g + geom_boxplot(aes(fill = Class)) +  stat_summary(fun.y = mean, geom = "line", 
    #                                                         lwd = 1, aes(group = Class, col = Class))
    # })
    
    output$allboxplot <- renderPlotly({
        df.m <- melt(data, id.vars = "Class")
        df.m <- df.m %>% drop_na()
        df.m$Class <- as.factor(df.m$Class)
        g <- ggplot(data = df.m, aes(x=variable, y=value)) 
        g <- g + geom_boxplot(aes(fill = Class)) +  stat_summary(fun.y = mean, geom = "line", 
                                                                 lwd = 1, aes(group = Class, col = Class))
        ggplotly(g)
    })
    
    output$corrplot <- renderPlotly({
        #GGally::ggcorr ((select(data, -Class)))
        corrtest <- cor(select(data, Clump_Thickness, Unif_Cell_Size, Unif_Cell_Shape, Marg_Adh, Single_Ep_Cell_Size, 
                               Bare_Nuclei, Bland_Chromatin, Normal_Nucleoli, Mitoses), method = "spearman")
        #corrplot::corrplot(corrtest)
        trace1 <- list(
            type = "heatmap", 
            x = colnames(select(data, -Class)), 
            y = colnames(select(data, -Class)), 
            z = corrtest
        )
        corrData <- list(trace1)
        layout <- list(title = "Features Correlation Matrix")
        p <- plot_ly()
        p <- add_trace(p, type=trace1$type, x=trace1$x, y=trace1$y, z=trace1$z)
        p <- layout(p, title=layout$title)
        p
    })
    
    # Downloadable csv of selected dataset ----
    output$downloadData <- downloadHandler(
        filename = function() {
            paste("BreastCancerDataset.csv", sep = "")
        },
        content = function(file) {
            write.csv(data, file, row.names = FALSE)
        }
    )
    
    output$cancerdata <- renderDataTable ({
        data
    })
    
    output$pcaresult <- renderDataTable ({
        pcRound <- round(PCs$rotation,3)
        datatable(as.data.frame(pcRound), rownames= TRUE)
    })
    
    output$biplot <- renderPlotly({
        simpleBiplot <- ggbiplot(PCs, choices = c(as.numeric(substr(input$pcs1,3,3)), 
                                                  as.numeric(substr(input$pcs2,3,3))), 
                                 ellipse=TRUE, labels=data$Class, groups = data$Class)
        ggplotly(simpleBiplot)
    })
    
    output$varpcplot <- renderPlot({
        par(mfrow = c(1, 2))
        plot(PCs$sdev^2/sum(PCs$sdev^2), xlab = "Principal Component",
             ylab = "Proportion of Variance Explained",  xlim = c(1,10), ylim = c(0, 1), type= 'b')
        plot(cumsum(PCs$sdev^2/sum(PCs$sdev^2)),xlab = "Principal Component",
             ylab = "Cum. Prop of Variance Explained", xlim = c(1,10), ylim = c(0, 1),type= 'b')
    })
    
    
    
    observe({
        
        validate(
            need(!is.null(input$checkGroup) , 
                 'Check at least one Predictor!')
        )
        
        set.seed(100)
        knn.pred <- knn(data.frame(trainData[,input$checkGroup]),
                        data.frame(testData[,input$checkGroup]),
                        trainData$Class, k = input$k)
        
        
        output$classError <- renderText({ paste("Classification Error = ",  round(ce(testData$Class, knn.pred), 3)) })
        
        output$confusionMatrix <- renderTable({
            # modify this to show title - confusion matrix
            # /false positive/positive false negative/negative
            true.positive    <- sum(knn.pred == 1 & testData$Class == 1)
            false.positive   <- sum(knn.pred == 1 & testData$Class == 0)
            true.negative    <- sum(knn.pred == 0 & testData$Class == 0)
            false.negative   <- sum(knn.pred == 0 & testData$Class == 1)
            row.names <- c("Prediction - FALSE", "Prediction - TRUE" )
            col.names <- c("Reference - FALSE", "Reference - TRUE")
            cbind(Outcome = row.names, as.data.frame(matrix(c(true.negative, false.negative, false.positive, true.positive) ,
                                                            nrow = 2, ncol = 2, dimnames = list(row.names, col.names))))
        }
        )
        
    })
    
    output$knnfrom10Fold <- renderPrint({
        knn.model
    })
    
    output$knnAccuracy <- renderPlotly ({
        f <- list(
            family = "Courier New, monospace",
            size = 18,
            color = "#7f7f7f"
        )
        x <- list(
            title = "K Neighbors",
            titlefont = f
        )
        y <- list(
            title = "Accuracy",
            titlefont = f
        )
        plot_ly(knn.model$results, x=knn.model$results$k, y=knn.model$results$Accuracy, mode = 'lines+markers') %>%
            layout(xaxis = x, yaxis = y)
    })
    
    output$classError10 <- renderText({ paste("Classification Error = ",  round(ce(testData$Class, knn.pred.cv), 3)) })
    
    output$confusionMatrix10 <- renderTable({
        # modify this to show title - confusion matrix
        # /false positive/positive false negative/negative
        truePositive    <- sum(knn.pred.cv == 1 & testData$Class == 1)
        falsePositive   <- sum(knn.pred.cv == 1 & testData$Class == 0)
        trueNegative    <- sum(knn.pred.cv == 0 & testData$Class == 0)
        falseNegative   <- sum(knn.pred.cv == 0 & testData$Class == 1)
        row.names <- c("Prediction - FALSE", "Prediction - TRUE" )
        col.names <- c("Reference - FALSE", "Reference - TRUE")
        cbind(Outcome = row.names, as.data.frame(matrix(c(trueNegative, falseNegative, falsePositive, truePositive) ,
                                                        nrow = 2, ncol = 2, dimnames = list(row.names, col.names))))
    })
    
    # observe({
    #     input_feature_x <- as.symbol(input$featureDisplay_x)
    #     input_feature_y <- as.symbol(input$featureDisplay_y)
    #     
    #     output$distPlotA <- renderPlotly ({
    #         # plot distribution of selected feature
    #         ggdistPlotA <- ggplot(data, aes_string(x = input$featureDisplay_x,
    #                                                fill = "factor(Class)")) +
    #             geom_histogram(position = "dodge")  +
    #             labs(x = input$featureDisplay_x,
    #                  y = "Count") + fte_theme() +
    #             scale_fill_manual(guide = F,values=c("#7A99AC", "#E4002B"))
    #         
    #     })
    #     
    #     output$distPlotB <- renderPlotly ({
    #         ggdistPlotB <- ggplot(data, aes_string(input$featureDisplay_y,
    #                                                fill = "factor(Class)")) +
    #             geom_histogram(position = "dodge") +
    #             labs(x = input$featureDisplay_y,
    #                  y = "Count") + fte_theme() +
    #             scale_fill_manual(guide = F,values=c("#7A99AC", "#E4002B"))
    #         
    #     })
    #     
    #     output$scatterPlotAB <- renderPlotly({
    #         # plot selected features against one another
    #         ggscatter <- ggplot(data, aes_string(x = input$featureDisplay_x, 
    #                                              y = input$featureDisplay_y, 
    #                                              color = "factor(Class)")) + 
    #             geom_point(size = 1, position = position_jitter(w = 0.1, h = 0.1)) + 
    #             labs(x = input$featureDisplay_x,
    #                  y = input$featureDisplay_y) +
    #             fte_theme() + 
    #             scale_color_manual(guide = F, values=c("#7A99AC", "#E4002B"))
    #         
    #     })
    #})
    
        output$knnAccPlot1 <- renderPlotly({
            f <- list(
                family = "Courier New, monospace",
                size = 18,
                color = "#7f7f7f"
            )
            x <- list(
                title = "K Neighbors",
                titlefont = f
            )
            y <- list(
                title = "Accuracy",
                titlefont = f
            )
            plot_ly(knn.model$results, x=knn.model$results$k, y=knn.model$results$Accuracy, mode = 'lines+markers') %>%
                layout(xaxis = x, yaxis = y)
        })
    
        
        v <- reactiveValues(doPlot = FALSE)
        
        observeEvent(input$plotAcc, {
            # 0 will be coerced to FALSE
            # 1+ will be coerced to TRUE
            v$doPlot <- input$plotAcc
        })
        
        observeEvent(input$checkGroup, {
            v$doPlot <- FALSE
        }) 
        
        output$knnAccPlot2 <- renderPlotly({
                
            if (v$doPlot == FALSE) return()
            
            isolate({
                
                trainSet <- data.frame(trainData[,input$checkGroup])
                testSet <-  data.frame(testData[,input$checkGroup])
                trainSet['Class'] <- trainData$Class
                
                set.seed(100)
                knn.model.upd <- train(Class ~., 
                                   data = trainSet, method = "knn",
                                   trControl=trctrl, 
                                   preProcess = c("center", "scale"), 
                                   tuneLength = 10, 
                                   tuneGrid = expand.grid(k = c(1:30)), 
                                   na.action = na.omit )
                
                f <- list(
                    family = "Courier New, monospace",
                    size = 18,
                    color = "#7f7f7f"
                )
                x <- list(
                    title = "K Neighbors",
                    titlefont = f
                )
                y <- list(
                    title = "Accuracy",
                    titlefont = f
                )
                plot_ly(knn.model.upd$results, x=knn.model.upd$results$k, 
                        y=knn.model.upd$results$Accuracy, mode = 'lines+markers') %>%
                    layout(xaxis = x, yaxis = y)
                
            })
        })
            
    
        observe({
            
            validate(
                need(!is.null(input$checkGroupLog) , 
                     'Check at least one Predictor!')
            )
            
            trainSet <- data.frame(trainData[,input$checkGroupLog])
            testSet <-  data.frame(testData[,input$checkGroupLog])
            
            glm.model <- glm(trainData$Class ~., data = trainSet, family = "binomial")
            glm.pred <- predict(glm.model, newdata = testSet, type = "response")
            glm.prob <- ifelse(glm.pred > 0.5, 1, 0)
            tbl1 <- table(testData$Class, glm.pred > 0.5)
            
            truePositive.log    <- sum(glm.prob == 1 & testData$Class == 1)
            falsePositive.log   <- sum(glm.prob == 1 & testData$Class == 0)
            trueNegative.log    <- sum(glm.prob == 0 & testData$Class == 0)
            falseNegative.log   <- sum(glm.prob == 0 & testData$Class == 1)
            
            
            output$logRegOut <- renderPrint({
                    summary(glm.model)
                
            })
            
            output$classErrorLog <- renderText({ str1 <- paste("Accuracy = ", round(sum(diag(tbl1))/sum(tbl1),3))
                                                 str2 <- paste("Classification Error = ",  
                                                               round(1 - sum(diag(tbl1))/sum(tbl1),3))
                                                 str3 <- paste("Sensitivity = ", round((truePositive.log/
                                                                                    (truePositive.log + falsePositive.log)),3) )
                                                 str4 <- paste("Specificity = ", round((trueNegative.log/
                                                                                    (trueNegative.log + falseNegative.log)),3))
                                                 paste(str1, str2, str3, str4, sep = "\n")
                                                 })
            
            output$confusionMatrixLog <- renderTable({
                # modify this to show title - confusion matrix
                # /false positive/positive false negative/negative
                
                row.names <- c("Prediction - FALSE", "Prediction - TRUE" )
                col.names <- c("Reference - FALSE", "Reference - TRUE")
                cbind(Outcome = row.names, as.data.frame(matrix(c(trueNegative.log, falseNegative.log, 
                                                                  falsePositive.log, truePositive.log) ,
                                                                nrow = 2, ncol = 2, dimnames = list(row.names, col.names))))
            })
        })
        
}

#shinyApp(ui, server)
