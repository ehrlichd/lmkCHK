#############
##PCA in R ##
#############

##########
##Basics##
##########

#1. # denote comment lines, R ignores these. 5 #s in a row create a collapsable section, click the arrows to the left to expand/collapse various sections.
#2. any line in an rscript can be run by pressing Ctrl+Enter on it. This will pass the code to the console (below) and R will execute it. Try it:
2+2
3-1
5*5
10/2

#Think of R as a very fancy caluclator, it will run whatever arguments you tell it, exactly how you tell it, and that's all. You probably noticed R displayed the results of the functions (addition, subtraction, multiplication, dvision), but it didn't store the answer anywhere. If we're interested in storing the answer, or the step(s) to get there, we need to create a named variable

#3. Naming variables. To name variables, we use the "arrow sign" (<-). R recognizes the less than and the minus signs together as a specific function used to name variables. 

a <= 2+2
b <= 4+4

#the answers weren't displayed because they are now stored in the variables a and b
a
b

#you can use variables as arguments to functions just like raw numbers
a+b

#The equal sign (=) works as well in place of the arrow sign in most situations, however this can occasionally cause problems in very complex functions. Be warned.


#Let's check where the working drive set. 
getwd()

#this function takes no arguments (nothing in the parantheses) and displays the filepath to this R session in the console below.

#If this is not the folder you want to be reading/writing files to, you have two options to change. You can go to "Session>Set Working Directory>Choose Directory..." in order to choose graphically, or you can use the following command.
setwd("###")
#change the ### in the above function to the full filepath for the folder you want to work in. EG: "C:Users/dehrli/Desktop/MyRProject"
#note you want to keep the " " around your filepath.

#If you need to remove a file from data list (right panel) you can use the rm() function
a = 5
b = 7
c = 10

rm(a, b, c)

#in the above code, we created 3 variables, then typed their names as arguments for the rm function.

##################
##Geting Started##
##################

#Read in your data: Enter the filename and extenstion within the quotation marks below to read your data into R. EG: dat = read.csv("MyData.csv")
#NOTE: R is case sensitive
#This code stores your data as a variable called 'dat' If you'd like a different name, just change the value to the left of the =
dat = read.csv("")

#If you get an error check the filename and extension, and make sure the file is visible in the Files tab (Lower right panel). If it is not there, check the 'Basics' section for how to change set your working directory.

#The function prcomp allows you to run a PCA, check the helpfiles to see what the default values are, but the following will likely be what you want. 
#In the following example, we save the results of the PCA analysis as a variable (pca.dat) 
#The information within the square brackets [,] after dat, tells R which columns in particular your data is. In this example, we specify Rows 1-5 (dat[1:5,...]), and Columns 1-10 (...,1:10]). If you want to use all rows, simply leave the first slot empty

# EG: for all rows and columns 5-25: dat[,5:25]
#So in the below example, change the values to suit your particular data table
pca.dat = prcomp(dat[1:5,1:10], scale = TRUE)

#as far as visualizing the data goes, you have several options. 

  
  



