#Greys,YlGnBu,Greens,YlOrRd,Bluered,RdBu,Reds,Blues,Picnic,Rainbow,Portland,Jet,Hot,Blackbody,Earth,Electric,Viridis,Cividis

library(plotly)
fes=read.table("fes.dat")
p<-plot_ly(x=fes[,1],y=fes[,2],z=fes[,3]/4.184,type="contour",autocontour = T,
  contours = list(
    start = 0, end = 15, size = 1), contours = list(showlabels = TRUE))%>%
  colorbar(title = "kcal/mol",titlefont=list(size=20,family="Droid Sans"))%>%
  layout(title = "",titlefont=list(size=30,family="Droid Sans"),
         xaxis = list(tickmode="auto",nticks=10,title = "CV-1 (nm)",tickfont=list(size=20),titlefont=list(size=30,family="Droid Sans"),range=c(0.8,1.5)),
         yaxis = list(tickmode="auto",nticks=10,title = "CV-2 (rad)",tickfont=list(size=20),titlefont=list(size=30,family="Droid Sans"),range=c(-3.14,3.14)))
p		

