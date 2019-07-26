library(caret)
library(RColorBrewer)

#Read Data and massage
data <- read_delim("C:/Users/Avisek/Downloads/ST558/Project3/breast-cancer-wisconsin.data", delim = ",", col_names = FALSE)
colnames(data) <- c("id_number", "Clump_Thickness", "Unif_Cell_Size", "Unif_Cell_Shape", "Marg_Adh", "Single_Ep_Cell_Size", "Bare_Nuclei", "Bland_Chromatin", "Normal_Nucleoli", "Mitoses", "Class")
data[data == '?'] <- NA
data <- data %>% drop_na()
data$Bare_Nuclei <- as.numeric(data$Bare_Nuclei)
data <- data %>% select(-id_number)
data$Class <- ifelse(data$Class == 4, 1, 0)
data$Class <- as.factor(data$Class)

#PCA Analysis
  PCs <- prcomp(select(data, -Class), center = TRUE, scale = TRUE)
  
  
  #Split Train and Test
  set.seed(100)
  trainIndex <- sample(1:nrow(data), size = nrow(data) * 0.8 )
  testIndex<- setdiff(1:nrow(data), trainIndex)
  
  trainData <- data[trainIndex, ]
  testData <- data[testIndex, ]

  
  #Train the kNN Model
  set.seed(100)
  trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
  set.seed(100)
  knn.model <- train(Class ~., 
                     data = trainData, method = "knn",
                     trControl=trctrl, 
                     preProcess = c("center", "scale"), 
                     tuneLength = 10, 
                     tuneGrid = expand.grid(k = c(1:30)), 
                     na.action = na.omit )
  set.seed(100)
  
  #Predict using test set
  knn.pred.cv <- predict(knn.model, 
                         newdata = dplyr::select(testData, -Class))
  
  
  
  fte_theme <- function() {
    
    # Generate the colors for the chart procedurally with RColorBrewer
    palette <- brewer.pal("Greys", n=9)
    color.background = palette[2]
    color.grid.major = palette[3]
    color.axis.text = palette[9]
    color.axis.title = palette[9]
    color.title = palette[9]
    text.size <- 14
    
    # Begin construction of chart
    theme_bw(base_size=9) +
      
      # Set the entire chart region to a light gray color
      theme(panel.background=element_rect(fill=color.background, color=color.background)) +
      theme(plot.background=element_rect(fill=color.background, color=color.background)) +
      theme(panel.border=element_rect(color=color.background)) +
      
      # Format the grid
      theme(panel.grid.major=element_line(color=color.grid.major,size=.25)) +
      theme(panel.grid.minor=element_blank()) +
      theme(axis.ticks=element_blank()) +
      
      # Format the legend, but hide by default
      theme(legend.background = element_rect(fill=color.background)) +
      theme(legend.text = element_text(size=text.size,color=color.axis.title)) +
      theme(legend.title = element_text(size=text.size,color=color.axis.title)) +
      theme(legend.position = "bottom") +
      theme(legend.direction = "vertical") +
      # Set title and axis labels, and format these and tick marks
      theme(plot.title=element_text(color=color.title, size=text.size, vjust=1.25)) +
      theme(axis.text.x=element_text(size=text.size,color=color.axis.text)) +
      theme(axis.text.y=element_text(size=text.size,color=color.axis.text)) +
      theme(axis.title.x=element_text(size=text.size,color=color.axis.title, vjust=0)) +
      theme(axis.title.y=element_text(size=text.size,color=color.axis.title, vjust=1.25)) +
      
      # Plot margins
      theme(plot.margin = unit(c(0.35, 0.2, 0.3, 0.35), "cm"))
  }