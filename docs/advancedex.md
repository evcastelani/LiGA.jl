---
layout: page
title: Advanced Examples
navigation: 7
---

## Advanced Examples

Here we are considering an example of Liga in circle detection. 
This problem is an important problem with several applications in
Computational Vision. 

Roughly speaking, this problem consists in, given a binarized image defined by 
A={(x1,y1),...,(xn,yn)} points try to get a circle with center (a,b) and radious
r that contain p points of A. 

The more classical way to solve this problem is using Hough Transform. 
For a good review of this subject we recommend [1] and [2]. The ideia of Hough 
Transform for circle detection is use a voting array and discretizations of the 
desired parameters (center and radius). Some improvements depending on coordinates
changing. 

Using GA we can get an elegant way for solving this problem without discretizations
(for example, without discretization for radius r). Of course we want to show 
a more speed way to solve the problem but an elegant way. Even that Liga is not
optimized for specific calculations. So, be careful. 

In sequence we have the code. This code can be downloaded [here](GA_dect.tar.gz).

	#load basic packages
	using PyPlot, Liga
	layout(2)
 
	#Auxiliary function	
	function findind(tensor::Array{Float64,3})
	  	(m,n,p)=size(tensor)
	  	vlmax=maximum(tensor)
	  	for i=1:m,j=1:n,k=1:p
              if tensor[i,j,k]==vlmax
                  return i,j,k,vlmax
              end
		end
	end


	function main(filename)
		#reading the example file
		A=readdlm(filename)
		#setting the size of the image and defining main parameters
		m=300
    	n=300
		maxrad=150	
		npun=length(A[:,1])	
		Vot=zeros(m,n,maxrad)
		#embeding A in conformal space and stored in EA
    	EA=Array{cmultvec,2}(npun,1)
    	einf=cb(id, true, false, 1.0)
		ea12=cb(e12, false, false, 1.0)
    	eaid=cb(id, true, false, -1.0)
    	SS=cmultvec(1.0*(e1,false,true))
		raio=1.0
		centro=[1.0,1.0]
		for i=1:npun
			EA[i]=conformal([A[i,1],A[i,2]])
		end
    	#display(EA)
		#the main block of the idea
    	#for i=1:npun,j=i+1:npun,k=j+1:npun #conservative method
		for i=1:2:npun,j=i+1:3:npun,k=j+1:5:npun #speed up				
			X=EA[i]^EA[j]^EA[k]
			SS=dual(X)
			raio = sqrt(((SS⋅SS)/((SS⋅einf)⋅(SS⋅einf))).comp[1].scl) 
			centro = conftore(projection(SS, ea12) / (SS ⋅ eaid))
        	if 1.0<=centro[1]<=m && 1.0<=centro[2]<=n
				if 1.0<=abs(raio)<=maxrad
					centro=round.(Int,centro)
					raio=abs.(round(Int,raio))
					Vot[centro[1],centro[2],raio]=Vot[centro[1],centro[2],raio]+1
				end
			end		
		end		
		sol=zeros(3)
		valsol=0.0
		sol[1],sol[2],sol[3],valsol=findind(Vot)
		display(sol)
		display(valsol)
		draw_solution(A,sol[1],sol[2],sol[3])
	end

	function draw_solution(A,cfx,cfy,rf)
		angulo=[0:((2*pi)/360):2*pi+1;]
		x=zeros(length(angulo))
		y=zeros(length(angulo))
		for i=1:length(angulo)
	    	x[i]=(abs(rf)*cos(angulo[i])+cfx)
	   		y[i]=(abs(rf)*sin(angulo[i])+cfy)
		end
		plot(A[:,1],A[:,2],".")
		plot(x,y,color="red","-")
		ax=gca()
		ax[:axis]("equal")
	end

As result of example datafile-50-300-300-10.txt
we have

![figure_1.jpeg]

[1] Duda, Richard O., and Peter E. Hart. "Use of the Hough transformation to 
detect lines and curves in pictures." Communications of the ACM 15.1 (1972): 11-15.

[2] Illingworth, John, and Josef Kittler. "A survey of the Hough transform.
"Computer vision, graphics, and image processing 44.1 (1988): 87-116.
