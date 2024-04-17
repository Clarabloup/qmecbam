cat("\014") # Clear console
rm(list=ls())# Clean workspace

# Libraries
library(data.table)
library(dplyr)

# directories
dir.bam <- "C:/Users/cchabrillan/Documents/RQmec/Tests A1/R/Fichiers textes/"

# Discharge measurements 
Q.station <- 'Batiscan_d_Q.txt'

setwd(dir.bam)

#two gauge stations
wl.h1 <- 'Batiscan_d_h.txt'
wl.h2 <- 'Becancour_d_h.txt'

interpolation.input='cubic'
interpolation.model='linear' 
dt.stage.record=60     # in minutes
dt.model = 1           # in minutes

all.data=data=list(h1=data.frame(),h2=data.frame())
all.data$h1<- fread(file.path(wl.h1))
all.data$h2 <- fread(file.path(wl.h2))

# Initial date to capture the tide
date.start <- data.frame('year'=2009,
                         'month'="Aug",
                         'day'=15,
                         'hour'=0,
                         'minute'=0,
                         'second'=0)

# Final date to capture the tide
date.end <- data.frame('year'=2009,
                       'month'="Aug",
                       'day'=30,
                       'hour'=0,
                       'minute'=0,
                       'second'=0)

# Pre-treatment of input data
for(i in 1:length(all.data)){
  
  colnames(all.data[[i]]) <- c('Date','Heure','h')
  
  all.data[[i]] <- data.frame(all.data[[i]][,1:2],
                              stage=all.data[[i]][,3])
  
  # transform and add date in a different format
  all.data[[i]]$date=paste(all.data[[i]]$Date,all.data[[i]]$Heure)
  all.data[[i]]$date=as.POSIXct(all.data[[i]]$date, format="%d-%b-%Y %H:%M:%S")
  # Water level measurement during period a specific period of time
  data[[i]]<- all.data[[i]][dplyr::between(all.data[[i]]$date,
                                           as.POSIXct("2009-06-15 00:00:00"),
                                           as.POSIXct("2009-06-30 00:00:00")),]
  # Interpolation input if NA detected
  gap=any(as.numeric(diff(data[[i]]$date),units='mins')!=dt.stage.record)
  if(gap==TRUE){time.inter.input.data = seq(data[[i]]$date[21],
                                data[[i]]$date[nrow(data[[i]])],by=dt.stage.record/60)
    if(interpolation.input=='cubic'){
      h=spline(x=data[[i]]$date,y=data[[i]]$h,method = "fmm",xout=time.inter.input.data)
    }else if(interpolation.input=='linear'){
      h=approx(x=data[[i]]$date,y=data[[i]]$h,xout=time.inter.input.data)
    }
    # Rebuild data frame with input data interpolated 
    data[[i]] = data.frame(year=year(time.inter.input.data),
                           month=month(time.inter.input.data),
                           day=mday(time.inter.input.data),
                           hour=hour(time.inter.input.data),
                           minute=minute(time.inter.input.data),
                           second=second(time.inter.input.data),
                           h=h$y,
                           date=time.inter.input.data,
                           u.h=0)
  }else{
    
    # Add uncertainty at gauge station. Assumption negligible
    data[[i]]$u.h <- 0
  }
}

# Merge results from water level measurements
data.wl <- merge(data$h1,data$h2,by=c('year','month','day','hour','minute','second','date'),suffixes =c('1','2'))
data.wl <- data.wl %>% arrange(date)

# Model environment 

# Interpolation if time step of stage data is higher than model's
inter.required =any(dt.stage.record!=dt.model)
if(inter.required==T){
  # interpolation.model to get a time step defined at dt.model  
  time.inter.model =seq(data.wl$date[1],data.wl$date[nrow(data.wl)],by=dt.model*dt.stage.record)
  
  if(interpolation.model=='cubic'){
    h1=spline(x=data.wl$date,y=data.wl$h1,method = "fmm",xout=time.inter.model)
    h2=spline(x=data.wl$date,y=data.wl$h2,method = "fmm",xout=time.inter.model)
  }else if(interpolation.model=='linear'){
    h1=approx(x=data.wl$date,y=data.wl$h1,xout=time.inter.model)
    h2=approx(x=data.wl$date,y=data.wl$h2,xout=time.inter.model)
  }
  
  # Rebuild data frame with observations interpolated 
  data.model.inter = data.frame(year=year(time.inter.model),
                                month=month(time.inter.model),
                                day=mday(time.inter.model),
                                hour=hour(time.inter.model),
                                minute=minute(time.inter.model),
                                second=second(time.inter.model),
                                date=time.inter.model,
                                h1=h1$y,
                                u.h1=0,
                                h2=h2$y,
                                u.h2=0)
  print(data.model.inter)
}else{
  data.model.inter = data.wl
}

ADCP.Q.measurements <- data.frame(fread(paste0(Q.station)))
ADCP.Q.measurements$date=paste(ADCP.Q.measurements$Date,ADCP.Q.measurements$Heure)
ADCP.Q.measurements$date=as.POSIXct(ADCP.Q.measurements$date, format="%d-%b-%Y %H:%M:%S")
time.inter.model =seq(ADCP.Q.measurements$date[1],ADCP.Q.measurements$date[nrow(ADCP.Q.measurements)],by=dt.stage.record/60)
print (time.inter.model)
print(ADCP.Q.measurements)
# Merge water level and discharge measurements (if)
data.model.Q.wl <- merge(data.model.inter, 
                         ADCP.Q.measurements, 
                         by = c("date"),
                         all.x = T,
                         all.y = T)

# cal.data = data.model.Q.wl[which(!is.na(data.model.Q.wl$Q))[1]:
#                              rev(which(!is.na(data.model.Q.wl$Q)))[1],]
cal.data = data.model.Q.wl

# Replace NA by -9999 for BaM! (manage of missing values)
cal.data[is.na(cal.data)] <- -9999

#save calibration data
write.table(cal.data,
            file=paste0(dir.bam,'calibrationData.txt'),
            sep='\t ',
            quote=F,
            col.names = T,
            row.names = F)
