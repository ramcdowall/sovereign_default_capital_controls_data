#Sovereign Debt and Capital Controls
#By: Robert A. McDowall
#Empirical investigation of capital controls and default risk
###############################################################################
#Libraries
library(reshape)
library(dplyr)
library(plm)
library(pglm)
library(stargazer)
library(multiwayvcov)
library(lmtest)
library(tidyr)
library(reshape2)
library(data.table)
library(ggplot2)
library('latticeExtra')

#setwd()
###############################################################################
#Helper Functions
# turn off scientific notation except for big numbers
options(scipen = 9)
# set larger font size for qplot (default is 12)
theme_set(theme_gray(base_size = 18))

### functions for correct SEs in regression tables

# function to calculate corrected SEs for OLS regression 
cse = function(reg) {
  rob = sqrt(diag(vcovHC(reg, type = "HC1")))
  return(rob)
}

# clustered SEs, clustered on "group"... could also cluster on "time" 
# compute Stata-like degrees of freedom adjustment for number of groups
# See http://www.richard-bluhm.com/clustered-ses-in-r-and-stata-2/
clse = function(reg) { 
  # index(reg, "id") returns the id or entity variable vector 
  G = length(unique(index(reg,"id")))
  N = length(index(reg,"id"))
  dfa = (G/(G - 1))   # note Bluhm multiplies this by finite-sample df adjustment
  rob = sqrt(diag(dfa*vcovHC(reg, method="arellano", type = "HC1", 
                             cluster = "group")))
  return(rob)
}

# corrected SEs for IV regressions
ivse = function(reg) {
  rob = robust.se(reg)[,2]
  return(rob)
}
###########################
#Read data
Data_SOC = read.csv("data/mcdowall_sovdefcc_data_annual.csv")

#Pooled Regression
reg1b = plm(boi ~ Yield_mean  + Govt_GDP +  GDP_growth ,
            data = Data_SOC,
            #family = binomial('probit'),
            index = c("Country","Year"), model="pooling")
summary(reg1b)

#Country Fixed Effects Regression
reg2b = plm(boi ~ Yield_mean   + Govt_GDP+ GDP_growth,
            data = Data_SOC, 
            #family = binomial('probit'),
            index = c("Country","Year"), model="within")
summary(reg2b)

#Country  + Time Fixed Effects Regression
reg3b = plm(boi ~ Yield_mean   + Govt_GDP+ GDP_growth,
            data = Data_SOC, 
            #family = binomial('probit'),
            index = c("Country","Year"), model="within", effect = "twoways")
summary(reg3b)

#Pooled Regression
reg1c = plm(kai ~ Yield_mean +  Govt_GDP +  GDP_growth,
            data = Data_SOC,
            #family = binomial('probit'),
            index = c("Country","Year"), model="pooling")
summary(reg1c)

#Country Fixed Effects Regression
reg2c = plm(kai ~ Yield_mean + Govt_GDP + GDP_growth,
            data = Data_SOC, 
            #family = binomial('probit'),
            index = c("Country","Year"), model="within")
summary(reg2c)

#Country +Time Fixed Effects Regression
reg3c = plm(kai ~ Yield_mean   + Govt_GDP + GDP_growth,
            data = Data_SOC, 
            #family = binomial('probit'),
            index = c("Country","Year"), model="within", effect = "twoways")
summary(reg3c)


stargazer(reg1b, reg2b, reg3b, reg1c, reg2c, reg3c,
          se=list(clse(reg1b),clse(reg2b),clse(reg3b), clse(reg1c),clse(reg2c), clse(reg3c)), 
          title="Panel regressions, clustered SEs", type="text", 
          column.labels=c( "Pooled OLS", "Country FE", "Two Ways", "Pooled OLS", "Country FE", "Two Ways"), 
          df=FALSE, digits=4)



##########################################
###############################################################################
##Create Plot for inflow restrictions against bond yields

#Import Data
Data_Uribe = read.csv('data/cc_indidces_uribe_et_al.csv')
Data_Spread = read.csv('data/mcdowall_10yr_spread_monthly.csv')

#Create Monthly CC Indices
#Repeate dataframe 
Data_Uribe2 = as.data.frame(lapply(Data_Uribe, rep, 12))
#Create ID 
Data_Uribe2 =Data_Uribe2 %>% group_by(Country, Year) %>% mutate(id = row_number())
#Combine with year to create identical id to spread data
Data_Uribe2 = Data_Uribe2 %>% group_by(Country, Year, id) %>% mutate(Year_month = paste("X", Year, "M", id, sep =''))


#Reshape spread data for panel 
Data_Spread2 = reshape(Data_Spread,
                       
                       varying = c( "X1995M1",          "X1995M2"       ,   "X1995M3"  ,       
                                    "X1995M4"      ,    "X1995M5"        ,  "X1995M6",          "X1995M7"     ,     "X1995M8"  ,       
                                    "X1995M9"       ,   "X1995M10"      ,   "X1995M11",         "X1995M12"   ,      "X1996M1"   ,      
                                    "X1996M2"        ,  "X1996M3"      ,    "X1996M4"  ,        "X1996M5"   ,       "X1996M6"    ,     
                                    "X1996M7"          ,"X1996M8"      ,    "X1996M9"    ,      "X1996M10"  ,       "X1996M11"     ,   
                                    "X1996M12",         "X1997M1"    ,      "X1997M2"    ,      "X1997M3" ,         "X1997M4"      ,   
                                    "X1997M5"  ,        "X1997M6"   ,       "X1997M7"     ,     "X1997M8",          "X1997M9"       ,  
                                    "X1997M10"   ,      "X1997M11"  ,       "X1997M12"      ,   "X1998M1",          "X1998M2"         ,
                                    "X1998M3"     ,     "X1998M4"  ,        "X1998M5"        ,  "X1998M6"          ,"X1998M7"         ,
                                    "X1998M8"     ,     "X1998M9",          "X1998M10"       ,  "X1998M11"       ,  "X1998M12"        ,
                                    "X1999M1"      ,    "X1999M2"          ,"X1999M3"         , "X1999M4"       ,   "X1999M5"         ,
                                    "X1999M6"        ,  "X1999M7"          ,"X1999M8",          "X1999M9"       ,   "X1999M10",        
                                    "X1999M11"       ,  "X1999M12"       ,  "X2000M1",          "X2000M2"     ,     "X2000M3",         
                                    "X2000M4"          ,"X2000M5"        ,  "X2000M6"  ,        "X2000M7"     ,     "X2000M8"  ,       
                                    "X2000M9",          "X2000M10"     ,    "X2000M11" ,        "X2000M12"  ,       "X2001M1" ,        
                                    "X2001M2" ,         "X2001M3"     ,     "X2001M4"   ,       "X2001M5"  ,        "X2001M6"  ,       
                                    "X2001M7"  ,        "X2001M8"    ,      "X2001M9"    ,      "X2001M10",         "X2001M11"  ,      
                                    "X2001M12"  ,       "X2002M1"   ,       "X2002M2"     ,     "X2002M3"          ,"X2002M4"    ,     
                                    "X2002M5"    ,      "X2002M6"  ,        "X2002M7"      ,    "X2002M8"         , "X2002M9"     ,    
                                    "X2002M10"    ,     "X2002M11",         "X2002M12"      ,   "X2003M1"        ,  "X2003M2"      ,   
                                    "X2003M3"      ,    "X2003M4"         , "X2003M5"        ,  "X2003M6"       ,   "X2003M7"       ,  
                                    "X2003M8"        ,  "X2003M9"         , "X2003M10"         ,"X2003M11"      ,   "X2003M12"        ,
                                    "X2004M1"         , "X2004M2"        ,  "X2004M3",          "X2004M4"      ,    "X2004M5"         ,
                                    "X2004M6",          "X2004M7"       ,   "X2004M8" ,         "X2004M9"     ,     "X2004M10"        ,
                                    "X2004M11" ,        "X2004M12"    ,     "X2005M1" ,         "X2005M2"   ,       "X2005M3"         ,
                                    "X2005M4"  ,        "X2005M5"     ,     "X2005M6"   ,       "X2005M7"   ,       "X2005M8"         ,
                                    "X2005M9"   ,       "X2005M10"  ,       "X2005M11"  ,       "X2005M12",         "X2006M1"         ,
                                    "X2006M2"    ,      "X2006M3"  ,        "X2006M4"    ,      "X2006M5"         , "X2006M6",         
                                    "X2006M7"     ,     "X2006M8" ,         "X2006M9"     ,     "X2006M10"       ,  "X2006M11",        
                                    "X2006M12"     ,    "X2007M1",          "X2007M2"      ,    "X2007M3"       ,   "X2007M4"  ,       
                                    "X2007M5"       ,   "X2007M6"         , "X2007M7"       ,   "X2007M8"      ,    "X2007M9"   ,      
                                    "X2007M10"        , "X2007M11"        , "X2007M12"        , "X2008M1"      ,    "X2008M2"     ,    
                                    "X2008M3"         , "X2008M4"       ,   "X2008M5",          "X2008M6"    ,      "X2008M7"     ,    
                                    "X2008M8",          "X2008M9"       ,   "X2008M10" ,        "X2008M11"   ,      "X2008M12"      ,  
                                    "X2009M1" ,         "X2009M2"      ,    "X2009M3"   ,       "X2009M4"   ,       "X2009M5"        , 
                                    "X2009M6"  ,        "X2009M7"     ,     "X2009M8"    ,      "X2009M9"  ,        "X2009M10"        ,
                                    "X2009M11" ,        "X2009M12"  ,       "X2010M1"    ,      "X2010M2",          "X2010M3"         ,
                                    "X2010M4"    ,      "X2010M5"   ,       "X2010M6"      ,    "X2010M7"          ,"X2010M8"         ,
                                    "X2010M9"    ,      "X2010M10",         "X2010M11"     ,    "X2010M12"         ,"X2011M1"         ,
                                    "X2011M2"     ,     "X2011M3"          ,"X2011M4"       ,   "X2011M5"         , "X2011M6",         
                                    "X2011M7"      ,    "X2011M8"         , "X2011M9"        ,  "X2011M10"       ,  "X2011M11",        
                                    "X2011M12"      ,   "X2012M1"        ,  "X2012M2"         , "X2012M3"       ,   "X2012M4"  ,       
                                    "X2012M5"        ,  "X2012M6"       ,   "X2012M7",          "X2012M8"      ,    "X2012M9"   ,      
                                    "X2012M10"        , "X2012M11"     ,    "X2012M12",         "X2013M1"     ,     "X2013M2"    ,     
                                    "X2013M3",          "X2013M4"     ,     "X2013M5"  ,        "X2013M6"    ,      "X2013M7"     ,    
                                    "X2013M8" ,         "X2013M9"    ,      "X2013M10"  ,       "X2013M11"  ,       "X2013M12"     ,   
                                    "X2014M1"  ,        "X2014M2"   ,       "X2014M3"    ,      "X2014M4"  ,        "X2014M5"       ,  
                                    "X2014M6"   ,       "X2014M7"  ,        "X2014M8"     ,     "X2014M9" ,         "X2014M10"       , 
                                    "X2014M11"   ,      "X2014M12",         "X2015M1"      ,    "X2015M2",          "X2015M3"         ,
                                    "X2015M4"      ,    "X2015M5"          ,"X2015M6"        ,  "X2015M7",          "X2015M8",         
                                    "X2015M9"       ,   "X2015M10"        , "X2015M11"        , "X2015M12"),
                       v.names = "Yield",
                       timevar = "Year",
                       times = c("X1995M1",          "X1995M2"       ,   "X1995M3"  ,       
                                 "X1995M4"      ,    "X1995M5"        ,  "X1995M6",          "X1995M7"     ,     "X1995M8"  ,       
                                 "X1995M9"       ,   "X1995M10"      ,   "X1995M11",         "X1995M12"   ,      "X1996M1"   ,      
                                 "X1996M2"        ,  "X1996M3"      ,    "X1996M4"  ,        "X1996M5"   ,       "X1996M6"    ,     
                                 "X1996M7"          ,"X1996M8"      ,    "X1996M9"    ,      "X1996M10"  ,       "X1996M11"     ,   
                                 "X1996M12",         "X1997M1"    ,      "X1997M2"    ,      "X1997M3" ,         "X1997M4"      ,   
                                 "X1997M5"  ,        "X1997M6"   ,       "X1997M7"     ,     "X1997M8",          "X1997M9"       ,  
                                 "X1997M10"   ,      "X1997M11"  ,       "X1997M12"      ,   "X1998M1",          "X1998M2"         ,
                                 "X1998M3"     ,     "X1998M4"  ,        "X1998M5"        ,  "X1998M6"          ,"X1998M7"         ,
                                 "X1998M8"     ,     "X1998M9",          "X1998M10"       ,  "X1998M11"       ,  "X1998M12"        ,
                                 "X1999M1"      ,    "X1999M2"          ,"X1999M3"         , "X1999M4"       ,   "X1999M5"         ,
                                 "X1999M6"        ,  "X1999M7"          ,"X1999M8",          "X1999M9"       ,   "X1999M10",        
                                 "X1999M11"       ,  "X1999M12"       ,  "X2000M1",          "X2000M2"     ,     "X2000M3",         
                                 "X2000M4"          ,"X2000M5"        ,  "X2000M6"  ,        "X2000M7"     ,     "X2000M8"  ,       
                                 "X2000M9",          "X2000M10"     ,    "X2000M11" ,        "X2000M12"  ,       "X2001M1" ,        
                                 "X2001M2" ,         "X2001M3"     ,     "X2001M4"   ,       "X2001M5"  ,        "X2001M6"  ,       
                                 "X2001M7"  ,        "X2001M8"    ,      "X2001M9"    ,      "X2001M10",         "X2001M11"  ,      
                                 "X2001M12"  ,       "X2002M1"   ,       "X2002M2"     ,     "X2002M3"          ,"X2002M4"    ,     
                                 "X2002M5"    ,      "X2002M6"  ,        "X2002M7"      ,    "X2002M8"         , "X2002M9"     ,    
                                 "X2002M10"    ,     "X2002M11",         "X2002M12"      ,   "X2003M1"        ,  "X2003M2"      ,   
                                 "X2003M3"      ,    "X2003M4"         , "X2003M5"        ,  "X2003M6"       ,   "X2003M7"       ,  
                                 "X2003M8"        ,  "X2003M9"         , "X2003M10"         ,"X2003M11"      ,   "X2003M12"        ,
                                 "X2004M1"         , "X2004M2"        ,  "X2004M3",          "X2004M4"      ,    "X2004M5"         ,
                                 "X2004M6",          "X2004M7"       ,   "X2004M8" ,         "X2004M9"     ,     "X2004M10"        ,
                                 "X2004M11" ,        "X2004M12"    ,     "X2005M1" ,         "X2005M2"   ,       "X2005M3"         ,
                                 "X2005M4"  ,        "X2005M5"     ,     "X2005M6"   ,       "X2005M7"   ,       "X2005M8"         ,
                                 "X2005M9"   ,       "X2005M10"  ,       "X2005M11"  ,       "X2005M12",         "X2006M1"         ,
                                 "X2006M2"    ,      "X2006M3"  ,        "X2006M4"    ,      "X2006M5"         , "X2006M6",         
                                 "X2006M7"     ,     "X2006M8" ,         "X2006M9"     ,     "X2006M10"       ,  "X2006M11",        
                                 "X2006M12"     ,    "X2007M1",          "X2007M2"      ,    "X2007M3"       ,   "X2007M4"  ,       
                                 "X2007M5"       ,   "X2007M6"         , "X2007M7"       ,   "X2007M8"      ,    "X2007M9"   ,      
                                 "X2007M10"        , "X2007M11"        , "X2007M12"        , "X2008M1"      ,    "X2008M2"     ,    
                                 "X2008M3"         , "X2008M4"       ,   "X2008M5",          "X2008M6"    ,      "X2008M7"     ,    
                                 "X2008M8",          "X2008M9"       ,   "X2008M10" ,        "X2008M11"   ,      "X2008M12"      ,  
                                 "X2009M1" ,         "X2009M2"      ,    "X2009M3"   ,       "X2009M4"   ,       "X2009M5"        , 
                                 "X2009M6"  ,        "X2009M7"     ,     "X2009M8"    ,      "X2009M9"  ,        "X2009M10"        ,
                                 "X2009M11" ,        "X2009M12"  ,       "X2010M1"    ,      "X2010M2",          "X2010M3"         ,
                                 "X2010M4"    ,      "X2010M5"   ,       "X2010M6"      ,    "X2010M7"          ,"X2010M8"         ,
                                 "X2010M9"    ,      "X2010M10",         "X2010M11"     ,    "X2010M12"         ,"X2011M1"         ,
                                 "X2011M2"     ,     "X2011M3"          ,"X2011M4"       ,   "X2011M5"         , "X2011M6",         
                                 "X2011M7"      ,    "X2011M8"         , "X2011M9"        ,  "X2011M10"       ,  "X2011M11",        
                                 "X2011M12"      ,   "X2012M1"        ,  "X2012M2"         , "X2012M3"       ,   "X2012M4"  ,       
                                 "X2012M5"        ,  "X2012M6"       ,   "X2012M7",          "X2012M8"      ,    "X2012M9"   ,      
                                 "X2012M10"        , "X2012M11"     ,    "X2012M12",         "X2013M1"     ,     "X2013M2"    ,     
                                 "X2013M3",          "X2013M4"     ,     "X2013M5"  ,        "X2013M6"    ,      "X2013M7"     ,    
                                 "X2013M8" ,         "X2013M9"    ,      "X2013M10"  ,       "X2013M11"  ,       "X2013M12"     ,   
                                 "X2014M1"  ,        "X2014M2"   ,       "X2014M3"    ,      "X2014M4"  ,        "X2014M5"       ,  
                                 "X2014M6"   ,       "X2014M7"  ,        "X2014M8"     ,     "X2014M9" ,         "X2014M10"       , 
                                 "X2014M11"   ,      "X2014M12",         "X2015M1"      ,    "X2015M2",          "X2015M3"         ,
                                 "X2015M4"      ,    "X2015M5"          ,"X2015M6"        ,  "X2015M7",          "X2015M8",         
                                 "X2015M9"       ,   "X2015M10"        , "X2015M11"        , "X2015M12"),
                       direction = "long")

#Keep only spread data
Data_Spread3 = Data_Spread2 %>% filter(Bond == "10YEARSPREAD")
Data = left_join(Data_Uribe2, Data_Spread3, by = c("Country" = "Country.Name", "Year_month" = "Year") )

#Filter out large open economies 
Data_SOC = Data %>% filter ( ! ( Country == "Germany" | Country == "Italy"|
                                   Country == "France"| Country == "United Kingdom" | Country == "United States" | Country == "Japan"
                                 | Country == "World Bank" | Country == "" ))

Data_SOC = Data_SOC %>% filter (!(incomegroup=="High income"))


Data_SOC$bondinflow = ifelse(Data_SOC$boi == 0, 0, 1 )
Data_SOC$inflow = ifelse(Data_SOC$kai == 0, 0, 1 )


###
mean_ts = Data_SOC %>% group_by(Year_month) %>% 
  summarise(mean_binfl = mean(boi, na.rm = TRUE), mean_infl = mean(kai, na.rm = TRUE), mean_yield = mean(Yield, na.rm = TRUE) )
mean_asts = mean_ts
mean_asts$Date <- gsub("X", "", mean_asts$Year_month)
mean_asts$Date <- gsub("M", "-", mean_asts$Date)
mean_asts$Date = as.character(mean_asts$Date, format='%Y-%m')
mean_asts$Date <- paste(mean_asts$Date, "-01", sep="")
mean_asts$Date = as.Date(mean_asts$Date, format = "%Y-%m-%d")

mean_asts$mean_binfl <- gsub("NaN", NA, mean_asts$mean_binfl)

mean_asts = within(mean_asts, rm(Year_month))



mean_asts$mean_binfl <- as.numeric(as.character(mean_asts$mean_binfl))


mean_asts = mean_asts %>% arrange(by = Date)


obj1 <- xyplot(mean_yield ~ Date, mean_asts, type = "l" , lwd=2, xlab = "Date",
               ylab = "Mean 10-Year Spread")
obj2 <- xyplot(mean_infl ~ Date, mean_asts, type = "l", lwd=3, xlab = "Date",
               ylab = "Mean Inflow Restrictions Index")
obj3 <- xyplot(mean_binfl ~ Date, mean_asts, type = "l", lwd=2)

# --> Make the plot with second y axis:
#setEPS()
#postscript("obj.eps")
doubleYScale(obj1, obj2, add.ylab2 = TRUE)
#saveas(gca, 'obj.eps','epsc');
dev.off()

