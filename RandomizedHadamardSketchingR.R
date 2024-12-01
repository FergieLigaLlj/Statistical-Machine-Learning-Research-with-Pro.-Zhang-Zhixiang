library(phangorn)
# Read the CSV file into a data frame
data <- read.csv("wind.csv",header=FALSE)

# Convert the data frame to a matrix
matrix_data <- as.matrix(data)

Y = matrix_data[,1]
X = matrix_data[,2:15]
p <- ncol(X)
### if the sample size is not to the power of two, padding X,y with zeros makes them comformable with SRHT
padding<-function(X,y){
  m<-nrow(X)
  if(ceiling(log(m,2))>log(m,2)){
    m1<-floor(log(m,2))+1
    padX<-rbind(X,matrix(0,2^m1-m,p))
    pady<-append(y,rep(0,2^m1-m))
  }
  else{padX<-X;pady<-y}
  return(list(padX=padX,pady=pady))
}


# Input: sketch size m, coefficient vector c, data X, y
# Output: sketched least square solutions r1, linear combination of r1, sketched data S*X and S*y
Esticoef_SRHT<-function(m,c,X,y,partial=0){
  p<-ncol(X)
  pad<-padding(X,y)
  X1<-pad$padX;y1<-pad$pady
  n1<-nrow(X1)
  gamma<-m/n1
  SXy<-apply(sample(c(1,-1),n1,replace=TRUE,prob=c(0.5,0.5))*cbind(X1,y1),2,fhm)[which(rbinom(n1,1,gamma)!=0),]/sqrt(m)
  if(partial==0){g<-qr.solve(SXy[,1:p],SXy[,p+1])}
  else{SX<-SXy[,1:p];M<-solve(t(SX)%*%SX);g<-M%*%(t(X)%*%y)}
  return(list(r1=g,r2=sum(c*g),r3=SXy))
}




start_time <- Sys.time()
m = 2000
c = 1:14
Esticoef_SRHT(m,c,X,Y,partial=0)
end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)