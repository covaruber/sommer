nonlinearGraphFixture <- function(){
  nodes <- list(
    x=list(op="input",column="x",center=1,scale=2),
    z=list(op="input",column="z"),
    hidden=list(op="affine",inputs=c("x","z"),
      weights=matrix(c(1,2,1,3),nrow=2),bias=c(4,4)),
    activation=list(op="tanh",input="hidden"),
    hidden2=list(op="affine",inputs="activation",
      weights=matrix(c(5,6,7,8),nrow=2),bias=c(9,10)),
    activation2=list(op="tanh",input="hidden2"),
    output=list(op="affine",inputs="activation2",
      weights=matrix(c(11,12),ncol=1),bias=13),
    identity=list(op="identity",input="output"))
  nonlinear_graph(nodes,output="identity",
    parameters=c("shared","rate","xweight","hiddenBias","w5","w6","w7",
      "w8","bias2","bias3","out1","out2","outBias"),
    parameterMap=list(shared=list(type="fixed",transform="log"),
      rate=list(type="fixed"),xweight=list(type="fixed"),
      hiddenBias=list(type="fixed"),w5=list(type="fixed"),w6=list(type="fixed"),
      w7=list(type="fixed"),w8=list(type="fixed"),bias2=list(type="fixed"),
      bias3=list(type="fixed"),out1=list(type="fixed"),out2=list(type="fixed"),
      outBias=list(type="fixed")))
}

test_that("native graph gradients and directional contractions match differences",{
  graph <- nonlinearGraphFixture()
  data <- data.frame(x=c(-1,0,1,2),z=c(2,1,0,-1),group=factor(c(1,1,2,2)))
  theta <- c(shared=log(.7),rate=.3,xweight=-.2,hiddenBias=.1,w5=.4,w6=-.3,
    w7=.2,w8=.5,bias2=-.1,bias3=.2,out1=.8,out2=-.4,outBias=-.15)
  direction <- matrix(seq(-.2,.19,length.out=3*length(theta)),
    nrow=length(theta),ncol=3)
  value <- nonlinear_graph_eval(graph,data,theta,direction)
  expect_length(value$output,nrow(data))
  expect_equal(dim(value$gradient),c(nrow(data),length(theta)))
  expect_equal(value$directionalFirst[,1],
    as.numeric(value$gradient%*%direction[,1]),tolerance=1e-12)
  expect_equal(value$directionalSecond[,1,2],value$directionalSecond[,2,1],
    tolerance=1e-12)
  step <- 1e-5
  for(j in seq_along(theta)){
    upper <- lower <- theta
    upper[j] <- upper[j]+step
    lower[j] <- lower[j]-step
    numericGradient <- (nonlinear_graph_eval(graph,data,upper)$output-
      nonlinear_graph_eval(graph,data,lower)$output)/(2*step)
    expect_equal(value$gradient[,j],numericGradient,tolerance=2e-7)
  }
  upper <- lower <- theta
  upper <- upper+step*direction[,1]
  lower <- lower-step*direction[,1]
  upperGradient <- nonlinear_graph_eval(graph,data,upper)$gradient
  lowerGradient <- nonlinear_graph_eval(graph,data,lower)$gradient
  expect_equal(value$hessianVector,
    (upperGradient-lowerGradient)/(2*step),tolerance=2e-6)
  upper <- lower <- theta
  upper <- upper+step*direction[,2]
  lower <- lower-step*direction[,2]
  upperHvp <- nonlinear_graph_eval(graph,data,upper,
    cbind(direction[,1],direction[,2],direction[,3]))$hessianVector
  lowerHvp <- nonlinear_graph_eval(graph,data,lower,
    cbind(direction[,1],direction[,2],direction[,3]))$hessianVector
  expect_equal(value$thirdDirectionalGradient,
    (upperHvp-lowerHvp)/(2*step),tolerance=2e-5)
  upper <- lower <- theta
  upper <- upper+step*direction[,3]
  lower <- lower-step*direction[,3]
  upperSecond <- nonlinear_graph_eval(graph,data,upper,direction)$directionalSecond
  lowerSecond <- nonlinear_graph_eval(graph,data,lower,direction)$directionalSecond
  expect_equal(value$directionalThird[,1,2,3],
    (upperSecond[,1,2]-lowerSecond[,1,2])/(2*step),tolerance=2e-5)
  expect_identical(unserialize(serialize(graph,NULL)),graph)
  expect_equal(nonlinear_graph_eval(graph,data,theta)$output,
    nonlinear_graph_eval(unserialize(serialize(graph,NULL)),data,theta)$output)
  saved <- tempfile(fileext=".rds")
  on.exit(unlink(saved),add=TRUE)
  saveRDS(graph,saved)
  script <- sprintf("graph <- readRDS(%s); stopifnot(inherits(graph, 'sommer_nonlinear_graph'), identical(graph$version, 1L))",
    encodeString(saved,quote='"'))
  status <- system2(file.path(R.home("bin"),"Rscript"),
    c("-e",shQuote(script)),stdout=TRUE,stderr=TRUE)
  expect_null(attr(status,"status"))
})

test_that("identity graph and tanh composition preserve declared semantics",{
  linear <- nonlinear_graph(list(x=list(op="input",column="x"),
    y=list(op="affine",inputs="x",weights=matrix(1L),bias=2L),
    output=list(op="identity",input="y")),output="output",
    parameters=c("slope","intercept"))
  data <- data.frame(x=c(-2,0,3))
  expect_equal(nonlinear_graph_eval(linear,data,
    c(slope=.4,intercept=1.2))$output,1.2+.4*data$x)
  logistic <- nonlinear_graph(list(
    x=list(op="input",column="x"),
    one=list(op="input",column="zero",center=-1),
    linear=list(op="affine",inputs="x",weights=matrix(1L),bias=2L),
    activation=list(op="tanh",input="linear"),
    mean=list(op="affine",inputs=c("activation","one"),
      weights=matrix(c(3L,3L),ncol=1),bias=0),
    output=list(op="identity",input="mean")),output="output",
    parameters=c("halfRate","shift","halfAsymptote"))
  x <- seq(-4,5,length.out=17)
  amplitude <- 3.2
  rate <- 1.4
  midpoint <- .6
  graphParameters <- c(halfRate=rate/2,shift=-rate*midpoint/2,
    halfAsymptote=amplitude/2)
  evaluated <- nonlinear_graph_eval(logistic,data.frame(x=x,zero=0),
    graphParameters)$output
  classic <- sommer_nonlinear_curve_eval("logistic",x,
    c(log(amplitude),log(rate),midpoint),gradient=FALSE)$mean
  expect_equal(evaluated,as.numeric(classic),tolerance=1e-12)
})

test_that("graph validation, group maps and memory guards reject invalid requests",{
  graph <- nonlinear_graph(list(x=list(op="input",column="x"),
    y=list(op="affine",inputs="x",weights=matrix(1L),bias=0L),
    output=list(op="identity",input="y")),output="output",
    parameters="u",parameterMap=list(u=list(type="random",group="group")))
  data <- data.frame(x=1:4,group=factor(c("a","a","b","b")))
  values <- matrix(c(.1,.1,.2,.2),ncol=1,
    dimnames=list(NULL,"u"))
  expect_length(nonlinear_graph_eval(graph,data,values)$output,4L)
  groupDirection <- array(0,c(4,1,3))
  groupDirection[1:2,1,1] <- 1
  directional <- nonlinear_graph_eval(graph,data,values,groupDirection)
  expect_equal(directional$directionalFirst[,1],c(1,2,0,0))
  groupDirection[2,1,1] <- .5
  expect_error(nonlinear_graph_eval(graph,data,values,groupDirection),
    "constant within group")
  values[2,1] <- .3
  expect_error(nonlinear_graph_eval(graph,data,values),"constant within group")
  expect_error(nonlinear_graph(list(x=list(op="relu",column="x")),
    output="x",parameters="p"),"Unsupported")
  expect_error(nonlinear_graph(list(x=list(op="identity",input="later"),
    later=list(op="input",column="x")),output="x",parameters="p"),
    "earlier graph node")
  expect_error(nonlinear_graph(list(x=list(op="input",column="x"),
    y=list(op="affine",inputs="x",weights=matrix(2L),bias=0L)),
    output="y",parameters="p"),"parameter indices")
  modified <- graph
  modified$version <- 99L
  expect_error(nonlinear_graph_eval(modified,data,values),"versioned object")
  expect_error(nonlinear_graph_eval(graph,data,values,
    control=sommer:::.mmes_control_internal(execution=list(denseMemoryMB=.00001))),"workspace")
})
