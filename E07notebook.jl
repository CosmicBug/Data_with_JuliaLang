### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# ╔═╡ b1608d70-5a3a-11eb-3caa-f9979f3c7a7b
begin
    using Pkg
    using LsqFit
    using CSV
    using PlutoUI
    using Plots
    using LsqFit
    import Pkg; Pkg.add("ORCA")
    using ORCA
    plotly()
    theme(:bright)
end

# ╔═╡ 3c6cb5d0-5a6b-11eb-2b6d-f5abd35ca0ee
begin
	using DataFrames
end

# ╔═╡ 1a31cd92-5b68-11eb-1690-ff8409553426
begin
	using Statistics
end

# ╔═╡ f63c9230-5b53-11eb-1238-93c94398f936
begin
	md"# E7e Magnetic Fields in Coils
	Shashank Shetty Kalavara #3728130
	"
end

# ╔═╡ aed75df0-5b60-11eb-3482-ed1fe9b038e5
md"Exporting data"

# ╔═╡ a9e44f30-5a69-11eb-331c-5737be2e247f
begin
	csv_file=CSV.File("E7_F1.dat")
	data_collective = DataFrame(csv_file)
	csv_file2=CSV.File("E7_L3.dat")
	data_collective2 = DataFrame(csv_file2)
	csv_file3=CSV.File("E7_P1_2R.dat")
	data_collective3 = DataFrame(csv_file3)
	csv_file4=CSV.File("E7_P1_R.dat")
	data_collective4 = DataFrame(csv_file4)
	csv_file5=CSV.File("E7_P1_R2.dat")
	data_collective5 = DataFrame(csv_file5)
end

# ╔═╡ 78877ab0-5a6a-11eb-3f62-791f69bc41c3
begin
	z=data_collective[:,1]./1000
	b=data_collective[:,2].*1.03 ./1000
	z2=data_collective2[:,1]./1000
	b2=data_collective2[:,2].*1.03 ./1000
	z3=data_collective3[:,1]./1000
	b3=data_collective3[:,2].*1.03 ./1000
	z4=data_collective4[:,1]./1000
	b4=data_collective4[:,2].*1.03 ./1000
	z5=data_collective5[:,1]./1000
	b5=data_collective5[:,2].*1.03 ./1000
end

# ╔═╡ a8f01762-5b6a-11eb-0575-7719b586fbe6


# ╔═╡ 74a82db0-5b67-11eb-26bd-e90a21eeb606
md"###### All error calculation are done through statistical standard deviation of the data "

# ╔═╡ 39ef3830-5b67-11eb-0fb7-a108615469a6
Total_Uncertainity = (0.00103+0.001)^(0.5) #through sqrt of least possible value measureable 

# ╔═╡ a7be8720-5b68-11eb-01c2-a5efc750954d
Standard_Dev_Distance = std(z)

# ╔═╡ c44cb560-5b68-11eb-19ea-57d30a8a6bf3
Standard_Dev_Bfeild = std(b)

# ╔═╡ a3b08e10-5b6a-11eb-04db-c5f799d9463b


# ╔═╡ a4ab68d0-5b6a-11eb-1318-3711dafdf01b


# ╔═╡ 08cf1740-5aa7-11eb-39cd-63169d28f033
md"#### 1. Magnetic field distribution in a long solenoid"

# ╔═╡ c89d4e10-5aa8-11eb-38ff-63d99cf36653

begin
    μ0 = 4*pi*10^(-7)
	@. Bz(a, I,L, R,N) = μ0*0.5*I*(N/L)*(a/(R^2+a^2)^(1/2)+(L-a)/(R^2+(L-a)^2)^(1/2))
	B = Bz(z2.-0.283,4,0.16,0.013,150)
end

# ╔═╡ b8429c50-5b70-11eb-21ad-c744615e4bbb
md"Equation [8] from the assignment is used for plotting the calculated values in this section"

# ╔═╡ 141e1540-5b71-11eb-1c93-d731f5561236
md"$B_z=B\left(a\right)=\frac{\mu \:_0}{2}\cdot ln\left(\frac{a}{\sqrt{R^2+a^2}}+\frac{L-a}{\sqrt{R^2+\left(L-a\right)^2}}\right)$"

# ╔═╡ 8fe90320-5b61-11eb-1584-a702773a07b8
md"Measured & calculated magnetic induction of a long solenoid"

# ╔═╡ 33b06500-5a6f-11eb-1cdc-dd0eccb4d84d
begin 
	plot(z2,b2,label="Measured", title="Measured & calculated magnetic induction",ylabel="B field [T]",xlabel="Distace [m]" )
	plot!(z2,B,label="Calculated")
end

# ╔═╡ 10c01b00-5b62-11eb-2e35-173727023502
md"###### c. Determining the magnetic induction at the end of the solenoid and compare experimental and theoretical values."

# ╔═╡ 74d86240-5b63-11eb-16ed-c1df10221833
md"As the length of the solenoid is given [0.16 m]. We can use this fact to extract a good value from the given data. Two lines are drawn in the graph along the Y axis with a constact X axis value, the distance between these two lines is set to be $0.16 m$. Now another line is drawn with constant value in Y axis, to obsrve the difference in value of line intersection between the 'Measured' and 'Calculated' plots "

# ╔═╡ 3ad6f260-5b62-11eb-2fc0-738d49217d1b
begin 
	plot(z2,b2,label="Measured", title="Measured & calculated magnetic induction",ylabel="B field [T]",xlabel="Distace [m]",legend=:topleft )
	plot!(z2,B,label="Calculated")
	ref=0.283
	plot!(ones(length(b2))*ref,B,color="green",label="Const. 0.283 m")#gree lines
	plot!(ones(length(b2))*(ref+0.16),B,color="green",label="Const. 0.443 m")
	#difference in these two values is always => 0.16m
	plot!(z2,ones(length(z2))*0.0023,color="grey",label="Const. 0.0023 T")
	plot!(z2,ones(length(z2))*0.0017,color="black",label="Const. 0.0017 T")#black line
	
end

# ╔═╡ 24963b00-5b62-11eb-2c22-a17033dd42d9
md" Magnetic Induction at the end of the solenoid for the measured experimental value is $0.0017$"

# ╔═╡ 8a4324a0-5b66-11eb-1c57-7fd9a1c79ed8
md" Magnetic Induction at the end of the solenoid for the calculated value is $0.0023 [T]$"

# ╔═╡ 6a15c7f0-5aa8-11eb-0001-a13b6e9ec73b
md"Differene in there values is $0.0006 [T]$"

# ╔═╡ e3948ff0-5b69-11eb-07f4-4d02b1d9a696
md"As the above seen Magnetic induction values are constituted within the larget standard error range of $\pm 0.0019 [T]$ a difference of $0.0006$, though graphically prominent in the plot, is quite negligible."

# ╔═╡ d010d5b0-5b69-11eb-08a8-9788b0334967
std(B)

# ╔═╡ 9e4d5020-5b6a-11eb-1281-6788d7d242c6


# ╔═╡ a0974020-5b6a-11eb-1195-d9e19ca08e81


# ╔═╡ 59476ce0-5aac-11eb-2fd1-03679f986130
md"#### 2. Magnetic field distribution of a flat circular coil"

# ╔═╡ af30d560-5aac-11eb-0075-15b01e3741ea
begin
	@. Bz2(zz, I, R,N)=((μ0*I*N*R^2)/2)*(1/(zz^2 + R^2)^3)^(1/2)
	BB2=Bz2(z.-0.264,3,0.093,390)
end

# ╔═╡ 80381840-5b60-11eb-08fb-9187f0c73feb
md"###### b. Ploting both the measured and calculated magnetic induction"

# ╔═╡ 8eb34470-5b70-11eb-36b9-d5088bfef1fc
md"Equation [9] from the assignment is used for making the calculated values plot in this section"

# ╔═╡ a12a70b0-5b70-11eb-04e1-83a05a696698
md"$B\left(z\right)=\frac{\mu _0\cdot I\cdot N\cdot R^2}{2}\cdot \frac{1}{\sqrt{\left(z^2+R^2\right)^3}}$"

# ╔═╡ 3b62ea20-5a6f-11eb-316c-ed528add78db
begin
	plot(z,b,title="Measured & calculated magnetic induction",label="Measured",ylabel="B field [T]",xlabel="Distace [m]")# the circle solanoid/coil
	plot!(z,BB2,label="Calculated")
end

# ╔═╡ 562c4c20-5b5a-11eb-1940-ef5bcc3a5dd9
md" ###### c. Determining at which distance from the coil center the maximum magnetic induction has decreases to half its value and comparing experimental and theoretical values "

# ╔═╡ b6b69200-5b58-11eb-1a90-2749fd278f93
md"#### $B\left(z\right)=\frac{\mu _0\cdot I\cdot N\cdot R^2}{2}\cdot \frac{1}{\sqrt{\left(z^2+R^2\right)^3}}$"

# ╔═╡ e10c2dd0-5b5d-11eb-0852-9164d2d1b986
md"We can use Equation [6] from the assignment, and calculate the inverse with respect to the variable $z$ i.e the distance from the coil center. And equate $B=\frac{Max- B -value}{2}$"

# ╔═╡ 158c8310-5b55-11eb-34dd-753b11791e45
begin
	md"### $z=\sqrt{\sqrt[3]{\left(\frac{\mu \:_0\cdot \:N\cdot \:I\cdot \:R^2}{B_h}\right)^2}-R^2} $    $  $ $→$    $B\left(z\right)=\frac{B_h}{2}$"
end

# ╔═╡ f2325e22-5b5f-11eb-25ed-6b817634ca6a
md"Distace of $\frac{B_{max}}{2}$ from the experimental values with maximum B feild as "

# ╔═╡ 11892240-5b60-11eb-01e8-21f3ca176d97
maximum(b)

# ╔═╡ 1c949330-5b5c-11eb-1d79-4170e154f5d5
md"#### $z_{experimental}=0.07131 \pm 0.005 [m]$"

# ╔═╡ 45c812a0-5b60-11eb-2e33-31118f33d54d
md"Distace of $\frac{B_{max}}{2}$ from the calculated values with maximum B feild as "

# ╔═╡ 4ee9b820-5b60-11eb-050f-1d8e84c39ede
maximum(BB2)

# ╔═╡ 5512daee-5b5c-11eb-2e7b-b34081a547b8
md"#### $ z_{calculated} = 0.08392 \pm 0.005 [m]$"

# ╔═╡ ebecd160-5aad-11eb-1fa2-3906d2afa1c1
md"#### 3. Magnetic field distribution of a pair of flat circular coils"

# ╔═╡ 0acbbfa0-5b6d-11eb-1561-e999526f97c6
md"###### b. Ploting both the measured and calculated magnetic induction."

# ╔═╡ 2a15e760-5aae-11eb-1eb3-59521e8ebb3a
begin
	plot(z3,b3,title="Measured & calculated magnetic induction",ylabel="B field [T]",xlabel="Distace [m]",label="2.R")
	plot!(z4,b4,label="R")
	plot!(z5,b5,label="R/2")
end

# ╔═╡ 4ecf3800-5b70-11eb-29b5-3d09f41eb7cd
md"Equation [10] from the assignment is used to make the calculated value graphs in this section"

# ╔═╡ 0db03a90-5b70-11eb-2ad1-6d61c7a7412b
md"##### $B\left(z\right)=\frac{\mu \:_0\cdot \:I\cdot \:N\cdot \:R^2}{2}\cdot \:\left(\frac{1}{\sqrt[3]{R^2+\left(z-\frac{b}{2}\right)^2}}+\frac{1}{\sqrt[3]{R^2+\left(z+\frac{b}{2}\right)^2}}\right)$"

# ╔═╡ 57b76c20-5aae-11eb-02f3-657e72816d3c
begin
	
	@. B_pair(a,I,R,N,b) = 0.5*μ0*I*N*R^2*(1/(R^2+(a-(b/2))^2)^(3/2)+(1)/(R^2+(a+(b/2))^2)^(3/2))
end

# ╔═╡ bee03cd2-5ab1-11eb-2694-b3859b0cb12d
begin
	plot(z3,b3,label="Experimental values",title="Measured & calculated magnetic induction",ylabel="B field [T]",xlabel="Distace [m]",legend=:bottomright)
	plot!(z3,B_pair(z3.-0.36,3.,0.096,390,2*0.096),label="Calculated values")
end

# ╔═╡ f783fcc0-5b47-11eb-1012-5bf898cebf2b
begin
	plot(z5,b5,label="Experimental values",title="Measured & calculated magnetic induction",ylabel="B field [T]",xlabel="Distace [m]",legend=:bottomright)
	plot!(z5,B_pair(z5.-0.288,3.,0.096,390,0.5*0.096),label="Calculated values")
end

# ╔═╡ 929c1e8e-5b48-11eb-2039-7943a4a51bec
begin
	plot(z4,b4,title="Measured & calculated magnetic induction",ylabel="B field [T]",xlabel="Distace [m]",label="Measured values",legend=:bottomright)
	plot!(z4,B_pair(z4.-0.31,3.,0.096,390,0.096),label="Calculated values")
end

# ╔═╡ 7075bfe0-5b6d-11eb-0833-3101d7071cb9
md"###### c. Determine for the Helmholtz configuration (b = R) theoretically and experimentally the magnetic field range with a deviation of the magnetic induction of maximum 5% from its central value."

# ╔═╡ 21af02f0-5b49-11eb-3faa-eb31f62ebb22
begin
	maximum(b4) #0.009094900000000001
	maximum(B_pair(z4.-0.31,3.,0.096,390,0.096)) #0.010958710980423536
	plot(z4,b4,legend=:bottom,ylabel="B field [T]",xlabel="Distace [m]",title="Helmholtz configuration b=R",label="measured")
	plot!(z4,B_pair(z4.-0.31,3.,0.096,390,0.096),label="Calculated values")
	plot!(z4,(maximum(b4)-(0.05*maximum(b4)))*ones(length(z4)),color="black",label="0.0086401")
	plot!(z4,(maximum(B_pair(z4.-0.31,3.,0.096,390,0.096))-(0.05*maximum(B_pair(z4.-0.31,3.,0.096,390,0.096))))*ones(length(z4)),color="blue",label="0.0104107")
end

# ╔═╡ 4fe446f0-5b6f-11eb-0385-71efecf08ca4
md"###### magnetic field range of the 'measured values' 0.36 - 0.26 $= 0.1 \pm 0.01 [m]$ "

# ╔═╡ 608f3eb0-5b6f-11eb-0ec5-6571a6742874
md"###### magnetic field range of the 'Calculated values' 0.355 - 0.265 = $0.09 \pm 0.01 [m]$ "

# ╔═╡ Cell order:
# ╟─b1608d70-5a3a-11eb-3caa-f9979f3c7a7b
# ╟─f63c9230-5b53-11eb-1238-93c94398f936
# ╟─3c6cb5d0-5a6b-11eb-2b6d-f5abd35ca0ee
# ╟─1a31cd92-5b68-11eb-1690-ff8409553426
# ╟─aed75df0-5b60-11eb-3482-ed1fe9b038e5
# ╟─a9e44f30-5a69-11eb-331c-5737be2e247f
# ╟─78877ab0-5a6a-11eb-3f62-791f69bc41c3
# ╟─a8f01762-5b6a-11eb-0575-7719b586fbe6
# ╟─74a82db0-5b67-11eb-26bd-e90a21eeb606
# ╠═39ef3830-5b67-11eb-0fb7-a108615469a6
# ╠═a7be8720-5b68-11eb-01c2-a5efc750954d
# ╟─c44cb560-5b68-11eb-19ea-57d30a8a6bf3
# ╟─a3b08e10-5b6a-11eb-04db-c5f799d9463b
# ╟─a4ab68d0-5b6a-11eb-1318-3711dafdf01b
# ╟─08cf1740-5aa7-11eb-39cd-63169d28f033
# ╟─c89d4e10-5aa8-11eb-38ff-63d99cf36653
# ╟─b8429c50-5b70-11eb-21ad-c744615e4bbb
# ╟─141e1540-5b71-11eb-1c93-d731f5561236
# ╟─8fe90320-5b61-11eb-1584-a702773a07b8
# ╟─33b06500-5a6f-11eb-1cdc-dd0eccb4d84d
# ╟─10c01b00-5b62-11eb-2e35-173727023502
# ╟─74d86240-5b63-11eb-16ed-c1df10221833
# ╟─3ad6f260-5b62-11eb-2fc0-738d49217d1b
# ╟─24963b00-5b62-11eb-2c22-a17033dd42d9
# ╟─8a4324a0-5b66-11eb-1c57-7fd9a1c79ed8
# ╟─6a15c7f0-5aa8-11eb-0001-a13b6e9ec73b
# ╟─e3948ff0-5b69-11eb-07f4-4d02b1d9a696
# ╠═d010d5b0-5b69-11eb-08a8-9788b0334967
# ╟─9e4d5020-5b6a-11eb-1281-6788d7d242c6
# ╟─a0974020-5b6a-11eb-1195-d9e19ca08e81
# ╟─59476ce0-5aac-11eb-2fd1-03679f986130
# ╟─af30d560-5aac-11eb-0075-15b01e3741ea
# ╟─80381840-5b60-11eb-08fb-9187f0c73feb
# ╟─8eb34470-5b70-11eb-36b9-d5088bfef1fc
# ╟─a12a70b0-5b70-11eb-04e1-83a05a696698
# ╟─3b62ea20-5a6f-11eb-316c-ed528add78db
# ╟─562c4c20-5b5a-11eb-1940-ef5bcc3a5dd9
# ╟─b6b69200-5b58-11eb-1a90-2749fd278f93
# ╟─e10c2dd0-5b5d-11eb-0852-9164d2d1b986
# ╟─158c8310-5b55-11eb-34dd-753b11791e45
# ╟─f2325e22-5b5f-11eb-25ed-6b817634ca6a
# ╠═11892240-5b60-11eb-01e8-21f3ca176d97
# ╟─1c949330-5b5c-11eb-1d79-4170e154f5d5
# ╟─45c812a0-5b60-11eb-2e33-31118f33d54d
# ╟─4ee9b820-5b60-11eb-050f-1d8e84c39ede
# ╟─5512daee-5b5c-11eb-2e7b-b34081a547b8
# ╟─ebecd160-5aad-11eb-1fa2-3906d2afa1c1
# ╟─0acbbfa0-5b6d-11eb-1561-e999526f97c6
# ╟─2a15e760-5aae-11eb-1eb3-59521e8ebb3a
# ╟─4ecf3800-5b70-11eb-29b5-3d09f41eb7cd
# ╟─0db03a90-5b70-11eb-2ad1-6d61c7a7412b
# ╟─57b76c20-5aae-11eb-02f3-657e72816d3c
# ╟─bee03cd2-5ab1-11eb-2694-b3859b0cb12d
# ╟─f783fcc0-5b47-11eb-1012-5bf898cebf2b
# ╟─929c1e8e-5b48-11eb-2039-7943a4a51bec
# ╟─7075bfe0-5b6d-11eb-0833-3101d7071cb9
# ╟─21af02f0-5b49-11eb-3faa-eb31f62ebb22
# ╟─4fe446f0-5b6f-11eb-0385-71efecf08ca4
# ╟─608f3eb0-5b6f-11eb-0ec5-6571a6742874
