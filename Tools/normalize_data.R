#Creator: Alexa Kur. | Arnau Vich

#Year: 2017

##USAGE##

# Copy or import function to R #

normalize = function(data,samples.in.rows = F,to.abundance = T,transform = "log",move.log = "min"){
  if (samples.in.rows == F){
    data = t(data)
    data[is.na(data)] = 0
  } else {
    data= as.matrix(data)
  }
  
  if(to.abundance == TRUE){
    data = sweep(data,1,rowSums(data),"/")
    data[is.nan(data)] = 0
  }
  
  
  if(!is.element(transform,c("log","asinsqrt","none"))) stop("Wrong transformation! Use 'log', 'asinsqrt' or 'none'")
  
  if (transform == "log"){
    if(move.log == "min"){
      min = min(data[data>0],na.rm = T)
      data = log(data+min)
      if(any(data<0,na.rm =T)) {data = data - min(data)}
    } else if(is.numeric(move.log)){
      data = log(data + move.log)
      if(any(data<0,na.rm =T)) {data = data - min(data)}
      
    } else{
      stop("please define the value to prevent -Inf from log transformation. can be numeric or 'min'")
    }
  } else if (transform == "asinsqrt" & to.abundance ==T) {
    data =asin(sqrt(data))
  } else if (transform == "asinsqrt" & any(data>1)) {
    stop("using asinsqrt without abundance transformation is not.allowed")
  } else {
    data = data
  }
  return(t(data))
} 
