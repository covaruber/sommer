multinomialModeFixture <- function(){
  list(
    response=c(1, 2, 3, 1, 3),
    predictor=rbind(c(-.3,.5),c(.2,-.4),c(1,.7),c(-.6,.1),c(.3,-.8)),
    design=Matrix::Matrix(rbind(c(1,0,.2),c(1,0,-.3),c(0,1,.5),c(0,1,.1),
      c(1,1,-.2),c(.1,1,0),c(.4,1,0),c(-.2,0,1),c(.2,0,1),c(.1,.3,1)),
      sparse=TRUE),
    precision=Matrix::Matrix(matrix(c(2,.2,-.1,.2,1.4,.15,-.1,.15,1.7),3),sparse=TRUE))
}

multinomialModeContext <- function(fixture, memory=1000){
  sommer:::sommer_multinomial_mode_create(fixture$response, fixture$design,
    fixture$precision,ncol(fixture$predictor),memory)
}
