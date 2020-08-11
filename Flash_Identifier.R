#!/usr/bin/env Rscript
if(!require(methods)){install.packages("methods")}
if(!require(argparser)){install.packages("argparser")}
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(zoo)){install.packages("zoo")}
if(!require(cli)){install.packages("cli")}
if(!require(log4r)){install.packages("log4r")}
p <- arg_parser("The flash calculator")
# Add a positional argument
p <- add_argument(p, "input_path", help="name of the file to analyze")
# Add a debug flag
p <- add_argument(p, "--debug", help="enable debug mode", flag=TRUE)
# Add a hill filter coefficient flag
p <- add_argument(p, "--filter", help="set a peak filter coefficient", default=0.3)
# Add a environmental background
p <- add_argument(p, "--en", help="set an environmental background", default=0)
# Add a environmental background
p <- add_argument(p, "--en_median", help="set an environmental background with median value of the valleys", flag=TRUE)
# Add a environmental background
p <- add_argument(p, "--en_minimum", help="set an environmental background with minimum value of the valleys", flag=TRUE)
argv <- parse_args(p)

if ((isTRUE(argv$en_median)) & (isTRUE(argv$en_minimum))){
stop("Please choose one among --en, --en_median and --en_minimum")
}

# function to display debug output
debug_msg <- function(...){
  if (isTRUE(argv$debug)){
    cat(paste0("DEBUG: ",...,"\n"))
  }
}

### set log file
logger <- create.logger()
logfile(logger) <- "flash_ide.log.txt"
level(logger) <- "DEBUG"

### create tryerror function
tryerror = function(x){
  tryCatch({
    sink("/dev/null")
    x
    sink()},
    error = function(msg) {
      error(logger,msg)
    },
    warning = function(msg) {
      warn(logger,msg)
    }
  ) 
  
}

###create needed functions
#valley
valley <- function(Fa,hill)
{valley <- c()
hill <- c(0,hill,length(Fa))
for (i in 2:length(hill))
{valley_candidate <- which(Fa[hill[i-1]:hill[i]]==min(Fa[hill[i-1]:hill[i]]))
valley[i-1] <- valley_candidate[length(valley_candidate)]+(hill[i-1]-1)
}
return (valley)
}
#hill
hill <- function(x){ 
  L <- length(x) 
  which( ((x[1:(L-2)]<x[2:(L-1)]) & (x[2:(L-1)]>=x[3:L]) ) 
  ) + 1 
} 
#subsequentneighbor
subsequentneighbor <- function(x,y){
  L <- y-x
  b <- which(L<0)
  y[b]=Inf
  return(y[which.min(y)])
}
#previousneighbor
previousneighbor <- function(x,y){
  L <- y-x
  b <- which(L>0)
  y[b]=-Inf
  return(y[which.max(y)])
}
#stralignR
stralignR <- function(x){
L <- c()
align <- c()
for (i in 1:length(x)){
L[i]  <- nchar(x[i])
}
for (i in 1:length(x)){
align[i] <- paste0(c(rep(" ",L[which.max(L)[1]]-L[i]),x[i]),collapse = "")
}
return (align)
}

###Read File
debug_msg("Debug option is set")
debug(logger,sprintf("Your input file is :%s ",argv$input_path))
debug_msg("Your input file is : ", argv$input_path)
tryerror(dirname(argv$input_path))
csv_dir_path = dirname(argv$input_path)
tryerror(sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(argv$input_path)))
csv_file_name = sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(argv$input_path))
tryerror(file.path(csv_dir_path, paste(csv_file_name, '_', '_result.txt', sep="")))
output_file_path = file.path(csv_dir_path, paste(csv_file_name, '_', '_result.txt', sep=""))

tryerror(read.table(argv$input_path, header=TRUE, sep=","))
Data = read.table(argv$input_path, header=TRUE, sep=",")
debug(logger,"Loading input file...")
debug_msg("Loading input file...")
tryerror(na.omit(Data))
Data = na.omit(Data)
debug(logger,sprintf("Omiting %.0f missing values...",sum(is.na(Data))))
debug_msg(sprintf("Omiting %.0f missing values...",sum(is.na(Data))))
tryerror(as.vector(Data[,1]))
Time <- as.vector(Data[,1])
debug(logger,"Loading time series...")
debug_msg("Loading time series...")
tryerror(as.vector(Data[,2]))
Fa <- as.vector(Data[,2])
debug(logger,"Loading value series...")
debug_msg("Loading value series...")
tryerror(data.frame(Fa,Time))
clean_data <- data.frame(Fa,Time)
tryerror(ggplot(clean_data,aes(Time,Fa))+geom_line()+xlab("Time")+ylab("Value")+ggtitle(paste(csv_file_name,"_original")))
overview <- ggplot(clean_data,aes(Time,Fa))+geom_line()+xlab("Time")+ylab("Value")+ggtitle(paste(csv_file_name,"_original"))
tryerror(ggsave(file.path(csv_dir_path, paste(csv_file_name,'_original.png'))))
debug(logger,"Saving original plot...")
debug_msg("Saving original plot...")



tryerror(hill(Fa))
Large <- hill(Fa)
debug(logger,"Identifying hills...")
debug_msg("Identifying hills...")
tryerror((quantile(Fa)[4]-quantile(Fa)[2])*argv$filter+quantile(Fa)[2])  ##²M°£fake peak
Large_threshold <- (quantile(Fa)[4]-quantile(Fa)[2])*argv$filter+quantile(Fa)[2]
debug(logger,sprintf("Setting hill filter:%f...",argv$filter))
debug_msg(sprintf("Setting hill filter:%f...",argv$filter))
tryerror(ggplot()+geom_line(aes(x=Time,y=Fa),color="grey")+geom_point(aes(x=Time[Large],y=Fa[Large]),shape=21, color="black", fill="#69b3a2")+ggtitle(paste(csv_file_name ,"_peaknotrim") )+xlab("Time")+ylab("Value"))
hillnotrim <- ggplot()+geom_line(aes(x=Time,y=Fa),color="grey")+geom_point(aes(x=Time[Large],y=Fa[Large]),shape=21, color="black", fill="#69b3a2")+ggtitle(paste(csv_file_name ,"_peaknotrim") )+xlab("Time")+ylab("Value")
tryerror(ggsave(file.path(csv_dir_path, paste(csv_file_name,'_peaknotrim.png'))))
debug(logger,"Saving peaknotrim plot...")
debug_msg("Saving peaknotrim plot...")



tryerror(seq(0.1,2,0.1))
coefficient <- seq(0.1,2,0.1)
tryerror(c())
hillcount <- c()
for (i in 1:length(coefficient)) {
  tryerror((quantile(Fa)[4]-quantile(Fa)[2])*coefficient[i]+quantile(Fa)[2])
  Large_threshold_test <- (quantile(Fa)[4]-quantile(Fa)[2])*coefficient[i]+quantile(Fa)[2] 
  tryerror(Large[which(Fa[Large]>Large_threshold_test)])
  Large_test <- Large[which(Fa[Large]>Large_threshold_test)]
  tryerror(length(Large_test))
  hillcount[i] <- length(Large_test)
}
tryerror(hillcount_df <- data.frame(coefficient,hillcount))
tryerror(for (i in 1:nrow(hillcount_df)){
  hillcount_df$text[i] <- paste0(c("(",hillcount_df$coefficient[i],",",hillcount_df$hillcount[i],")"),collapse="")
})
tryerror(ggplot(hillcount_df,aes(x=coefficient,y=hillcount))+geom_line()+geom_point(shape=21, color="black", fill="#8CEA00")+ggtitle(paste(csv_file_name,"_coefficient-peakcount"))+xlab("Coefficient")+ylab("Number of peak remained")+geom_text(label=hillcount_df$text, nudge_x = 0.25, nudge_y = 0.25, check_overlap = T))
hillcount_plot <- ggplot(hillcount_df,aes(x=coefficient,y=hillcount))+geom_line()+geom_point(shape=21, color="black", fill="#8CEA00")+ggtitle(paste(csv_file_name,"_coefficient-peakcount"))+xlab("Coefficient")+ylab("Number of peak remained")+geom_text(label=hillcount_df$text, nudge_x = 0.25, nudge_y = 0.25, check_overlap = T)
tryerror(ggsave(file.path(csv_dir_path, paste(csv_file_name,'_coefficient-peakcount.png'))))
debug(logger,"Saving coefficient-peakcount plot...")
debug_msg("Saving coefficient-peakcount plot...")

tryerror(Large[which(Fa[Large]>Large_threshold)])
Large <- Large[which(Fa[Large]>Large_threshold)]
tryerror(Time[Large][1:length(Time[Large])-1])
Time_Large_0 = Time[Large][1:length(Time[Large])-1]
tryerror(Time[Large][2:length(Time[Large])])
Time_Large_1 = Time[Large][2:length(Time[Large])]
tryerror(Time_Large_1-Time_Large_0)
Time_Large_index = Time_Large_1-Time_Large_0

if (length(which(Time_Large_index<mean(Time_Large_index)*0.25)) != 0)
{tryerror(which(Time_Large_index<mean(Time_Large_index)*0.25))
cluster_index = which(Time_Large_index<mean(Time_Large_index)*0.25)
remove_index = c()
for (i in c(1:length(cluster_index)))
{tryerror(c(Fa[Large[cluster_index[i]]],Fa[Large[cluster_index[i]+1]]))
smaller_index = c(Fa[Large[cluster_index[i]]],Fa[Large[cluster_index[i]+1]])
tryerror(which.min(smaller_index))
smaller = which.min(smaller_index)
if (smaller == 1)
{
tryerror(cluster_index[i])
remove_index[i] <- cluster_index[i]
}else
{
tryerror(cluster_index[i]+1)
remove_index[i] <- cluster_index[i]+1}
}
tryerror(Large[-remove_index])
Large <- Large[-remove_index]
}else 
{Large = Large}
debug(logger,"Filtering peaks...")
debug_msg("Filtering peaks successfully")
tryerror(ggplot()+geom_line(aes(x=Time,y=Fa))+geom_point(aes(x=Time[Large],y=Fa[Large]),shape=21, color="black", fill="#69b3a2")+geom_hline(yintercept = Large_threshold,col='#424200' )+ggtitle(paste(csv_file_name ,"_peaktrimed"))+xlab("Time")+ylab("Value"))
hilltrimed <- ggplot()+geom_line(aes(x=Time,y=Fa))+geom_point(aes(x=Time[Large],y=Fa[Large]),shape=21, color="black", fill="#69b3a2")+geom_hline(yintercept = Large_threshold,col='#424200' )+ggtitle(paste(csv_file_name ,"_peaktrimed"))+xlab("Time")+ylab("Value")
tryerror(ggsave(file.path(csv_dir_path, paste(csv_file_name,'_peaktrimed.png'))))
debug(logger,"Saving peaktrimed plot...")
debug_msg("Saving peaktrimed plot...")


tryerror(valley(Fa,Large))
Minimum <- valley(Fa,Large)
debug(logger,"Identifying valleys...")
debug_msg("Identifying valleys...")
tryerror(data.frame(start = c(0) , end = c(0)) )
flash_df <- data.frame(start = c(0) , end = c(0)) 
tryerror(Minimum[1])
flash_df[1,1] <- Minimum[1]
for (i in 1:length(Minimum))
{ 
  tryerror(subsequentneighbor(subsequentneighbor(flash_df[i,1],Large),Minimum))
  flash_df[i,2] <- subsequentneighbor(subsequentneighbor(flash_df[i,1],Large),Minimum)
  if (flash_df[i,2]==Minimum[length(Minimum)]) break
  tryerror(flash_df[i,2])
  flash_df[i+1,1] <- flash_df[i,2]
}
tryerror(ggplot()+geom_line(aes(x=Time,y=Fa))+geom_point(aes(x=Time[Large],y=Fa[Large]),shape=21, color="black", fill="#69b3a2")+geom_point(aes(x=Time[Minimum],y=Fa[Minimum]),shape=21, color="black", fill="#FFD306")+ggtitle(paste(csv_file_name,"_peak & valley")))
hillvalley_plot <- ggplot()+geom_line(aes(x=Time,y=Fa))+geom_point(aes(x=Time[Large],y=Fa[Large]),shape=21, color="black", fill="#69b3a2")+geom_point(aes(x=Time[Minimum],y=Fa[Minimum]),shape=21, color="black", fill="#FFD306")+ggtitle(paste(csv_file_name,"_peak & valley"))
tryerror(ggsave(file.path(csv_dir_path, paste(csv_file_name,'_peak & valley.png'))))
debug(logger,"Saving peak_valley plot...")
debug_msg("Saving peak & valley plot...")
tryerror(write.table(flash_df,file.path(csv_dir_path, paste(csv_file_name,'_flash table.csv'))))
debug(logger,"Saving flash table...")
debug_msg("Saving flash table...")


Fa_cluster <- list()
for (i in 1:(nrow(flash_df)))
{ tryerror(Fa[c(flash_df[i,1]):c(flash_df[i,2]) ])
  Fa_cluster[[i]] <- Fa[c(flash_df[i,1]):c(flash_df[i,2]) ] } 
Time_cluster <- list()
for (i in 1:(nrow(flash_df)))
{ tryerror(Time[c(flash_df[i,1]):c(flash_df[i,2]) ])
  Time_cluster[[i]] <- Time[c(flash_df[i,1]):c(flash_df[i,2]) ]}
debug(logger,"Building value clusters and time clusters...")
debug_msg("Building value clusters and time clusters...")


area_all <- c()
for (j in 1:length(Fa_cluster))
{area <- c()
for (i in 1:(length(Fa_cluster[[j]])-1))
{ tryerror((Fa_cluster[[j]][i+1]+Fa_cluster[[j]][i])*(Time_cluster[[j]][i+1]-Time_cluster[[j]][i])*0.5)
  area[i] <- (Fa_cluster[[j]][i+1]+Fa_cluster[[j]][i])*(Time_cluster[[j]][i+1]-Time_cluster[[j]][i])*0.5
}
tryerror(sum(area))
area_all[j] <- sum(area)
}
debug(logger,"Calculating curve area under each flash...")
debug_msg("Calculating curve area under each flash...")


tryerror(argv$en)
base_line = argv$en
if (isTRUE(argv$en_median)){
tryerror(median(Fa[Minimum]))
base_line = median(Fa[Minimum])
}else if (isTRUE(argv$en_minimum)){
tryerror(Fa[Minimum][which.min(Fa[Minimum])])
base_line = Fa[Minimum][which.min(Fa[Minimum])]  
}
debug(logger,sprintf("Enviromental background:%f",base_line))
debug_msg(sprintf("Enviromental background:%f",base_line))


base_area <- c()
for (j in 1:length(Time_cluster))
{ tryerror(base_line*(Time_cluster[[j]][length(Time_cluster[[j]])]-Time_cluster[[j]][1]))
  base_area[j] <- base_line*(Time_cluster[[j]][length(Time_cluster[[j]])]-Time_cluster[[j]][1])
}
tryerror(area_all-base_area)
flash_energy = area_all-base_area
tryerror(flash_energy)
flash_df$energy <- flash_energy
tryerror(write.table(flash_df,sep=",",row.names=F,file.path(csv_dir_path, paste(csv_file_name,'_flash table.csv'))))
debug_msg("Calculating flash energy...")
debug_msg("Saving flash table...")


peak <- c()
for (j in 1:length(Fa_cluster))
{ tryerror(max(Fa_cluster[[j]]))
  peak[j] <- max(Fa_cluster[[j]])  
}
tryerror(peak)
flash_df$peak <- peak
tryerror(write.table(flash_df,sep=",",row.names=F,file.path(csv_dir_path, paste(csv_file_name,'_flash table.csv'))))
debug(logger,"Saving flash table...")
debug_msg("Saving flash table...")


tryerror(seq(1,nrow(flash_df),1))
flash_df$number <- seq(1,nrow(flash_df),1)
tryerror(ggplot(flash_df,aes(x = number,y= peak)) +geom_bar(stat="identity", color="#e9ecef", alpha=0.9) +ggtitle(paste(csv_file_name ,"_peak")) +theme(plot.title = element_text(size=15) ))
peak_plot <- ggplot(flash_df,aes(x = number,y= peak)) +geom_bar(stat="identity", color="#e9ecef", alpha=0.9) +ggtitle(paste(csv_file_name ,"_peak")) +theme(plot.title = element_text(size=15) )
tryerror(ggsave(file.path(csv_dir_path, paste(csv_file_name,'_peak_barplot.png'))))
debug(logger,"Saving peak barplot...")
debug_msg("Saving peak barplot...")


tryerror(ggplot(flash_df,aes(x = number,y= energy)) +geom_bar(stat="identity", color="blue", alpha=0.9,fill="white") +ggtitle(paste(csv_file_name ,"_energy")) +theme(plot.title = element_text(size=15) ))
peak_plot <- ggplot(flash_df,aes(x = number,y= energy)) +geom_bar(stat="identity", color="blue", alpha=0.9,fill="white") +ggtitle(paste(csv_file_name ,"_energy")) +theme(plot.title = element_text(size=15) )
tryerror(ggsave(file.path(csv_dir_path, paste(csv_file_name,'_energy_barplot.png'))))
debug(logger,"Saving energy barplot...")
debug_msg("Saving energy barplot...")


tryerror(stralignR(c("max_peak",max(flash_df$peak))) )
max_peak <- stralignR(c("max_peak",max(flash_df$peak)))
tryerror(stralignR(c("median_peak",median(flash_df$peak))) )
median_peak <- stralignR(c("median_peak",median(flash_df$peak)))
tryerror(stralignR(c("max_energy",max(flash_df$energy))) )
max_energy <- stralignR(c("max_energy",max(flash_df$energy)))
tryerror(stralignR(c("median_energy",median(flash_df$energy))) )
median_energy <- stralignR(c("median_energy",median(flash_df$energy))) 
tryerror(stralignR(c("frequency",nrow(flash_df)/( Time[flash_df$end[nrow(flash_df)]]-Time[flash_df$start[1]]))) )
frequency <- stralignR(c("frequency",nrow(flash_df)/( Time[flash_df$end[nrow(flash_df)]]-Time[flash_df$start[1]])))
tryerror(cat("### File\n", argv$input_path, file=output_file_path, sep="\n", append=FALSE) )
tryerror(cat("\n### Descriptive statistics\n", file=output_file_path, append=TRUE) )

tryerror(c(max_peak,median_peak,max_energy,median_energy,frequency) )
statistics <-  c(max_peak,median_peak,max_energy,median_energy,frequency) 
tryerror(c(statistics[c(1,3,5,7,9)],"\n") )
stastistics_name <- c(statistics[c(1,3,5,7,9)],"\n") 
tryerror(statistics[c(2,4,6,8,10)] )
stastistics_value <- statistics[c(2,4,6,8,10)] 
tryerror(cat(stastistics_name,sep="   ",file=output_file_path, append=TRUE) )
tryerror(cat(stastistics_value,sep="   ",file=output_file_path, append=TRUE) )
debug(logger,"Outputing statistics...")
debug_msg("Outputing statistics...")


tryerror(cli_ol() )
tryerror(cli_li("File") )
tryerror(ulid <- cli_ul() )
tryerror(cli_li(argv$input_path) )
tryerror(cli_end(ulid) )
tryerror(cli_li("Statistics") )
tryerror(ulid <- cli_ul() )
tryerror(cli_li(paste0("Max_peak: ",as.character(max(flash_df$peak)))) )
tryerror(cli_li(paste0("Median_peak: ",as.character(median(flash_df$peak)))) )
tryerror(cli_li(paste0("Max_energy: ",as.character(max(flash_df$energy)))) )
tryerror(cli_li(paste0("Median_energy: ",as.character(median(flash_df$energy)))) )
tryerror(cli_li(paste0("Max_peak: ",as.character(nrow(flash_df)/( Time[flash_df$end[nrow(flash_df)]]-Time[flash_df$start[1]]      )))) )
tryerror(cli_end() )


debug(logger,"Finish")
message('Finish')





