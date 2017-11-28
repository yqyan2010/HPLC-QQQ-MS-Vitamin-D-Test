""" Replace Copy and Paste Raw VitD Data (in ng/mL) to Adjust """

""" Author: YYan"""
""" Date: Sep 19 2017 """

""" Load library packages """
import pandas as pd
from pandas import ExcelWriter
import numpy as np
import os
import re
#import xlsxwriter

""" Get lab staff name """
ana_by=input("Analyzed by: ")
rev_by=input("Reviewed by: ")
pub_by=input("Published by: ")

""" Get Sample Test Info """
test_day=str(input("Enter Test Date YYMMDD: "))

""" Get Sample Info """
num_pun3=int(input("Total Number of 3.2 Punch: "))
pun3={}
sp_pun3_name=[]
sp_pun3_num=[]
for i in range(0,num_pun3):
    id=str(input("Enter 3.2 Punch Sample ID: ")).upper()
    sp_pun3_name.append(id)
    num=int(input("Number of punches? "))
    sp_pun3_num.append(num)
    pun3[id]=num

""" Functions """
def find_date():
    import datetime
    date=str(datetime.date.today())
    return date

""" //Crate Library To Store Info// """ ## Library is Dictionary in Python
""" lib VSP target """
vsp_lib_d3={"VSP 93":61,"VSP 95":23,"VSP 97":11,"VSP 98":24,"VSP 103":33}
vsp_lib_d2={"VSP 93":np.nan,"VSP 97":17,"VSP 98":31,"VSP 103":3}

""" lib patient target """
pat_d3={"Patient 102":16,'Patient 90':32,'Patient 75':47,'Patient 66':25}
pat_d2={"Patient 102":29,'Patient 90':0,'Patient 75':0,'Patient 66':1}

""" lib ABVD """
ABVD_d3_lib={"ABVD-01":19.84,"ABVD-02":36.53,"ABVD-03":27.96,"ABVD-07":67.23,"ABVD-08":10.68,"ABVD-09":19.2}
ABVD_d2_lib={"ABVD-01":25.65,"ABVD-02":0.9,"ABVD-03":0.79,"ABVD-07":1.28,"ABVD-08":4.96,"ABVD-09":1.74}

""" lib punch size """
pun_fac={"6mm x1":1,"3.2mm x4":1.1378,"3.2mm x3":0.8533,"3.2mm x2":0.568889}

""" lib Calibration """
cal={"L1":1.953125,"L2":3.90625,"L3":7.8125,"L4":15.625,"L5":31.25,"L6":62.5,"L7":125,"L8":250}

""" //Create DataFrame//"""
""" Calibration, ABVD, Punch Data Frame"""
df_cal=(pd.DataFrame(data=cal,index=["exp"])).T
df_punch=(pd.DataFrame(data=pun_fac,index=["Factor"])).T
df_abvd=pd.DataFrame(data=ABVD_d3_lib,index=[""])
df_abvd=df_abvd.append(ABVD_d2_lib,ignore_index=True)
df_abvd.index=("25OHD3_Target","25OHD2_Target")
df_abvd=df_abvd.T
df_abvd["Total 25OHD"]=df_abvd["25OHD3_Target"]+df_abvd["25OHD2_Target"]

""" //Load Data// """
""" !!Imporant!! LOAD D2D3 ng/mL FROM LC/MS REPORTS """
# rpt parses first (data file name) and last (ng/mL) columns

#folder=r"D:\MassHunter\Data\Vitamin D"
#folder=folder+'\VITD_'+test_day+'\QuantReports\VITD_'+test_day

folder=r"C:\Users\yyan\Documents\Data\VitD_DBS\WKD"
folder=os.path.join(folder,'VITD_'+test_day)
f_path=os.path.join(folder,r"QuantReport_ISTD_ResultsSummary_FAST_B_04_00.xlsx")
rpt=pd.read_excel(f_path,sheetname="Summary",skiprows=7,names=["compound","ng/mL"],header=None,index_col=0,parse_cols="A,B,G")
lst=[i for i in range(0,len(rpt))]
rpt["Row"]=lst # Add a new column of numbers

name=rpt.loc['Target Compound']['compound'].iloc[0]
if name[5]=='3':
    d3_start=rpt.loc['Target Compound']['Row'].iloc[0]
    d2_start=rpt.loc['Target Compound']['Row'].iloc[1]
    rpt_d3=rpt[d3_start:d2_start-1]
    del rpt_d3['Row']
    del rpt_d3['compound']
    rpt_d2=rpt[d2_start:]
    del rpt_d2["Row"]
    del rpt_d2['compound']
else:
    d2_start=rpt.loc['Target Compound']['Row'].iloc[0]
    d3_start=rpt.loc['Target Compound']['Row'].iloc[1]
    rpt_d2=rpt[d2_start:d3_start-2]
    del rpt_d2['Row']
    del rpt_d2['compound']
    rpt_d3=rpt[d3_start:]
    del rpt_d3['Row']
    del rpt_d3['compound']

""" Build Info Table df_info """
#today=find_date() # find today's date
yr=test_day[0:2]
mon=test_day[2:4]
day=test_day[4:]
date=mon+"/"+day+"/"+"20"+yr
date_col=[date,date,date,date,date]
ext_id=yr+mon+day+"A"
ass_id=yr+mon+day+"B"
name_col=[ext_id,ass_id,ana_by,rev_by,pub_by]
ana_by_ini="".join(i for i in ana_by if i.isupper())
rev_by_ini="".join(i for i in rev_by if i.isupper())
pub_by_ini="".join(i for i in pub_by if i.isupper())
analyst_col=[ana_by_ini,ana_by_ini,ana_by_ini,rev_by_ini,pub_by_ini]
# Create df_info
df_info=pd.DataFrame(data=None,index=["Extracted ID","Assay ID","Data Analyzed By","Reviewed By","Published By"])
df_info["Name"]=name_col
df_info["Date"]=date_col
df_info["Analyst"]=analyst_col

""" //Quanlity Control (UTAK)//"""
""" UTAK D3 df_ut_d3 """
def build_df_ut_d3(a,b,c):
    raw_utL_d3=rpt_d3.loc[a].iloc[0] #"UTAK LW_1.d" analytes keep constant name
    raw_ut1_d3=rpt_d3.loc[b].iloc[0] #"UTAK L1_1.d" analytes keep constant name
    raw_ut2_d3=rpt_d3.loc[c].iloc[0] #"UTAK L2_1.d" keep constat name
    raw_ut_d3_col=[raw_utL_d3,raw_ut1_d3,raw_ut2_d3]
    # Build UTAK D3 Table
    df_ut_d3=pd.DataFrame(data=None,index=["Utak L","Utak 1","Utak 2"])
    df_ut_d3["D3 Expected"]=["9.3-12.7","23-32","63-85"]
    df_ut_d3["Target"]=[10,30,73]
    df_ut_d3["Raw (ng/mL)"]=raw_ut_d3_col
    df_ut_d3["Factor"]=list(df_ut_d3["Raw (ng/mL)"]/df_ut_d3["Target"])
    fac_avg_ut_d3=df_ut_d3["Factor"].mean()
    fac_std=df_ut_d3["Factor"].std() # fac_std" name will be re-used. Find standard deviation
    fac_rsd=fac_std/fac_avg_ut_d3 # Find Relative Standard Error
    df_ut_d3["Adjusted"]=[(df_ut_d3.loc[item,"Raw (ng/mL)"])/fac_avg_ut_d3 for item in df_ut_d3.index]
    df_ut_d3["Average"]=list([np.NaN,np.NaN,np.NaN])
    df_ut_d3.loc["Utak 2","Average"]=fac_avg_ut_d3 # Put fac_avg in Average Column and last row
    df_ut_d3["SD"]=list([np.NaN,np.NaN,np.NaN])
    df_ut_d3.loc["Utak 2","SD"]=fac_std # Put fac_std in "SD" column and last row
    df_ut_d3["CV or RSD"]=list([np.NaN,np.NaN,np.NaN])
    df_ut_d3.loc["Utak 2","CV or RSD"]=fac_rsd # Put fac_rsd in "CV or RSD" column and last row
    return df_ut_d3,fac_avg_ut_d3

""" UTAK D2 df_ut_d2 """
def build_df_ut_d2(a,b,c):
    raw_utL_d2=rpt_d2.loc[a].iloc[0] #"UTAK LW_1.d" keep constant name
    raw_ut1_d2=rpt_d2.loc[b].iloc[0] #"UTAK L1_1.d" keep constant name
    raw_ut2_d2=rpt_d2.loc[c].iloc[0] #"UTAK L2_1.d" keep constant name
    raw_ut_d2_col=[raw_utL_d2,raw_ut1_d2,raw_ut2_d2]
    # Build UTAK D2 Table
    df_ut_d2=pd.DataFrame(data=None,index=["Utak L","Utak 1","Utak 2"])
    df_ut_d2["D2 Expected"]=["8.5-11.5","25-34","62-84"]
    df_ut_d2["Target"]=[10,30,73]
    df_ut_d2["Raw (ng/mL)"]=raw_ut_d2_col
    df_ut_d2["Factor"]=list(df_ut_d2["Raw (ng/mL)"]/df_ut_d2["Target"])
    fac_avg_ut_d2=df_ut_d2["Factor"].mean()
    fac_std=df_ut_d2["Factor"].std() # Replace fac_std from D3 to D2
    fac_rsd=fac_std/fac_avg_ut_d2 # Replace fac_std from D3 to D2
    df_ut_d2["Adjusted"]=[df_ut_d2.loc[item,"Raw (ng/mL)"]/fac_avg_ut_d2 for item in df_ut_d2.index]
    df_ut_d2["Average"]=list([np.NaN,np.NaN,np.NaN])
    df_ut_d2.loc["Utak 2","Average"]=fac_avg_ut_d2
    df_ut_d2["SD"]=list([np.NaN,np.NaN,np.NaN])
    df_ut_d2.loc["Utak 2","SD"]=fac_std
    df_ut_d2["CV or RSD"]=list([np.NaN,np.NaN,np.NaN])
    df_ut_d2.loc["Utak 2","CV or RSD"]=fac_rsd
    return df_ut_d2,fac_avg_ut_d2

"""Create utak d3 dataframe"""
a='' # UTAK LW.d
b='' # UTAK L1.d
c='' # UTAK L2.d
rpt_d3_index=rpt_d3.index ## rpt_d3 and rpt_d2 may have different index
rpt_d2_index=rpt_d2.index
if 'UTAK LW.d' in rpt_d3_index:
    a='UTAK LW.d'
else:
    a='UTAK LW_1.d'
if 'UTAK L1.d' in rpt_d3_index:
    b='UTAK L1.d'
else:
    b='UTAK L1_1.d'
if 'UTAK L2.d' in rpt_d3_index:
    c='UTAK L2.d'
else:
    c='UTAK L2_1.d'
df_ut_d3,fac_avg_ut_d3=build_df_ut_d3(a,b,c)

"""Create utak d2 dataframe"""
if 'UTAK LW.d' in rpt_d2_index:
    a='UTAK LW.d'
else:
    a='UTAK LW_1.d'
if 'UTAK L1.d' in rpt_d2_index:
    b='UTAK L1.d'
else:
    b='UTAK L1_1.d'
if 'UTAK L2.d' in rpt_d2_index:
    c='UTAK L2.d'
else:
    c='UTAK L2_1.d'
df_ut_d2,fac_avg_ut_d2=build_df_ut_d2(a,b,c)

""" //Quanlity Control (VSP)// """
""" !!Important!! FIND CORRECTION FACTOR (Using VSP Data) """
""" Create VSP df_vsp Table"""
vsp_name=[]
for item in rpt_d3.index:
    if item.startswith("ADJ"):
        vsp_name.append(item)
vsp_name.sort()

df_vsp_d3=pd.DataFrame(data=None,index=[i.replace("ADJ","VSP")[:-2] for i in vsp_name])
df_vsp_d3["D3 ADJ Expected"]=np.nan

df_vsp_d2=pd.DataFrame(data=None,index=[i.replace("ADJ","VSP")[:-2] for i in vsp_name if "95" not in i if "93" not in i])
df_vsp_d2["D2 ADJ Expected"]=np.nan

df_vsp_d3["Target"]=[vsp_lib_d3[item[:-2]] for item in df_vsp_d3.index]
df_vsp_d2["Target"]=[vsp_lib_d2[item[:-2]] for item in df_vsp_d2.index]

df_vsp_d3["Raw"]=[rpt_d3.loc[item.replace("VSP","ADJ").__add__(".d")].iloc[0] for item in df_vsp_d3.index]
df_vsp_d2["Raw"]=[rpt_d2.loc[item.replace("VSP","ADJ").__add__(".d")].iloc[0] for item in df_vsp_d2.index]

fac_vsp_d3=[(df_vsp_d3.loc[item,"Raw"])/(df_vsp_d3.loc[item,"Target"]) for item in df_vsp_d3.index]
fac_vsp_d2=[(df_vsp_d2.loc[item,"Raw"])/(df_vsp_d2.loc[item,"Target"]) for item in df_vsp_d2.index]

fac_d3_avg=sum(fac_vsp_d3)/len(fac_vsp_d3) # ADJUST FACTOR (Average of factors)
fac_d2_avg=sum(fac_vsp_d2)/len(fac_vsp_d2) # ADJUST FACTOR (Average of factors)

df_vsp_d3["Adjusted"]=[(df_vsp_d3.loc[item,"Raw"])/(fac_d3_avg) for item in df_vsp_d3.index]
df_vsp_d2["Adjusted"]=[(df_vsp_d2.loc[item,"Raw"])/(fac_d2_avg) for item in df_vsp_d2.index]

df_vsp_d3["Factor"]=fac_vsp_d3 # Find Average value for adjust factors (adjustors)
df_vsp_d2["Factor"]=fac_vsp_d2

df_vsp_d3["Average"]=np.NaN # Store Average Factor D3 in a new Average column
row_num=len(df_vsp_d3)-1
col_num=list(df_vsp_d3.columns).index("Average")
df_vsp_d3.iloc[row_num,col_num]=fac_d3_avg ## Put average value in the last cell of Average Column
df_vsp_d2["Average"]=np.NaN # D2 Average Factor
row_num=len(df_vsp_d2)-1
col_num=list(df_vsp_d2.columns).index("Average")
df_vsp_d2.iloc[row_num,col_num]=fac_d2_avg

fac_d3_std=df_vsp_d3["Factor"].std() # Find D3 Standard Deviation of Adjust Factor
fac_d2_std=df_vsp_d2["Factor"].std() # D2 STD of Adjust Factor
df_vsp_d3["SD"]=np.NaN # D3 STD of Adjust factor in "SD" column last row
row_num=len(df_vsp_d3)-1
col_num=list(df_vsp_d3.columns).index("SD")
df_vsp_d3.iloc[row_num,col_num]=fac_d3_std
df_vsp_d2["SD"]=np.NaN # D2 STD of Adjust factor in "SD" column last row
row_num=len(df_vsp_d2)-1
col_num=list(df_vsp_d2.columns).index("SD")
df_vsp_d2.iloc[row_num,col_num]=fac_d3_std

df_vsp_d3["CV or RSD"]=np.NaN #"RSD=STD/AVG", D3 RSD in column and last row
row_num=len(df_vsp_d3)-1
col_num=list(df_vsp_d3.columns).index("CV or RSD")
df_vsp_d3.iloc[row_num,col_num]=(fac_d3_std)/(fac_d3_avg)
df_vsp_d2["CV or RSD"]=np.NaN # D2 RSD in column and last row
row_num=len(df_vsp_d2)-1
col_num=list(df_vsp_d2.columns).index("CV or RSD")
df_vsp_d2.iloc[row_num,col_num]=(fac_d2_std)/(fac_d2_avg)

del row_num,col_num

""" //Data Analysis// """
""" !!Important!! RAW DATA D2D3 ng/mL FROM LC-MS REPORT """
def find_sp_name(x):
    sp_name=[]
    if x.lower()=='d3':
        for item in rpt_d3_index: # rpt_d3 and rpt_d2 may have different index
            if item.startswith("B17"):
                sp_name.append(item)
            if item.startswith("Patient"):
                sp_name.append(item)
            if item.startswith("ABVD"):
                sp_name.append(item)
    sp_name.sort()
    if x.lower()=='d2':
        for item in rpt_d2_index: # rpt_d3 and rpt_d2 may have different index
            if item.startswith("B17"):
                sp_name.append(item)
            if item.startswith("Patient"):
                sp_name.append(item)
            if item.startswith("ABVD"):
                sp_name.append(item)
    sp_name.sort()
    return sp_name

# Raw D3 Values (ng/mL) from LC-MS report
sp_d3=[]
sp_name_d3=find_sp_name('d3')
for item in sp_name_d3:
    sp_d3.append(rpt_d3.loc[item].iloc[0])
# Raw D2 Values (ng/mL) from LC-MS report
sp_d2=[]
sp_name_d2=find_sp_name('d2')
for item in sp_name_d2:
    sp_d2.append(rpt_d2.loc[item].iloc[0])


""" //Data Analysis// """
""" !!Important!! D3D2 FINAL REPORTED VALUES HERE """
""" Create Sample Table df_sample and Calculation D3D2"""

if sp_name_d3 == sp_name_d2:
    df_sample=pd.DataFrame(data=None,index=[i[:-2] for i in sp_name_d3])
    df_sample["D3 Raw"]=sp_d3
    df_sample["D2 Raw"]=sp_d2
else:
    if len(sp_name_d3)> len(sp_name_d2):
        df_sample=pd.DataFrame(data=None,index=[i[:-2] for i in sp_name_d3])
        dif=list(set(sp_name_d3)-set(sp_name_d2))
        for item in dif:
            index=sp_name_d3.index(item)
            sp_d2.insert(index,0)
        df_sample["D3 Raw"]=sp_d3
        df_sample["D2 Raw"]=sp_d2
    if len(sp_name_d2)> len(sp_name_d3):
        df_sample=pd.DataFrame(data=None,index=[i[:-2] for i in sp_name_d2])
        dif=list(set(sp_name_d2)-set(sp_name_d3))
        for item in dif:
            index=sp_name_d2.index(item)
            sp_d3.insert(index,0)
        df_sample["D3 Raw"]=sp_d3
        df_sample["D2 Raw"]=sp_d2

""" 3.2 MM Punch Adjust Block (Optional) """
repeat=len(df_sample)
Pun_6=[1 for i in range(0,repeat)]
Pun_3=[0 for i in range(0,repeat)]

df_sample["# of Punchese 6.0mm"]=Pun_6
df_sample["# of Punchese 3.2mm"]=Pun_3

df_sample["Punch Adjusted D3 adj."]=df_sample["D3 Raw"]*(df_sample["# of Punchese 6.0mm"]+df_sample["# of Punchese 3.2mm"])
df_sample["Punch Adjusted D2 adj."]=df_sample["D2 Raw"]*(df_sample["# of Punchese 6.0mm"]+df_sample["# of Punchese 3.2mm"])

if num_pun3 != 0:
    df_sample_index=df_sample.index
    for item in df_sample_index:
        for i in range(0,num_pun3):
            if item.startswith(sp_pun3_name[i]):
                loc=list(df_sample_index).index(item) # Get location of 3.2mm adjusted sample in df_sample name index
                df_sample.loc[item,"# of Punchese 6.0mm"]=0
                df_sample.loc[item,"# of Punchese 3.2mm"]=sp_pun3_num[i]
                if sp_pun3_num[i]== 2:
                    df_sample.loc[item,"Punch Adjusted D3 adj."]=sp_d3[loc]/pun_fac["3.2mm x2"] # Update D3 adjusted value for 3.2mm spots
                    df_sample.loc[item,"Punch Adjusted D2 adj."]=sp_d2[loc]/pun_fac["3.2mm x2"]
else:
    pass

df_sample["D3 Final"]=[df_sample.loc[item,"Punch Adjusted D3 adj."]/fac_d3_avg/fac_avg_ut_d3 for item in df_sample.index] # D3 ng/ml FINAL REPORT VALUE
# Create D2 Final Column, remove D2 value <2
sp_d2_final=[df_sample.loc[item,"Punch Adjusted D2 adj."]/fac_d2_avg/fac_avg_ut_d2 for item in df_sample.index]
for i in range(0,len(sp_d2_final)):
    if sp_d2_final[i]<4:
        sp_d2_final[i]="<4" # Remove D2 final that <2
df_sample["D2 Final"]=sp_d2_final # D2 ng/mL FINAL REPORT Values
# Create 25OHD Final column, not include D2<2
OHD_final=[]
for item in df_sample.index:
    if not type(df_sample.loc[item,"D2 Final"]) == str:
        OHD_final.append(df_sample.loc[item,"D2 Final"]+df_sample.loc[item,"D3 Final"])
    else:
        OHD_final.append(df_sample.loc[item,"D3 Final"])
df_sample["Final 25OHD"]=OHD_final # TOTAL VitD ng/mL FINAL REPORT VAVLUE
# Create QC Recovery column
pat_rec=[]
rec=float()
for item in df_sample.index:
    if str(item).startswith('Patient'):
        if '_' in item:
            item_new=item[0:item.find('_')]
        else:
            item_new=item
        if item in pat_d3:
            rec=(df_sample.loc[item]['Final 25OHD']/(pat_d3[item_new]+pat_d2[item_new]))
            pat_rec.append(rec)
        else:
            pat_rec.append(np.NaN)
    else:
        pat_rec.append(np.NaN)
df_sample["Control Percent Recovery"]=pat_rec

""" //Quanlity Control (CAL)// """
""" Calibration Curve (VitD ng/mL) vs (VitD Concentration)"""
# Vit D3D2 Calibrated Value From LC-MS Repot
cal_name=[]
for i in range(1,9):
    cal_name.append("DS CAL L"+str(i)+".d") # "DS CAL L1-8.d" Constant Calibration Sample Name to Use
cal_d3=[]
for item in cal_name:
    if item in rpt_d3_index:
        cal_d3.append(rpt_d3.loc[item].iloc[0])
    else:
        cal_d3.append(np.NaN)
df_cal["25OHD3"]=cal_d3

cal_d2=[]
for item in cal_name:
    if item in rpt_d2_index:
        cal_d2.append(rpt_d2.loc[item].iloc[0])
    else:
        cal_d2.append(np.NaN)
df_cal["25OHD2"]=cal_d2

df_cal["Notes"]=np.NaN

""" //Write To Excel File// """
writer=pd.ExcelWriter("VitDresult.xlsx",engine='xlsxwriter')
df_info.to_excel(writer,"Sheet1",na_rep="")
df_ut_d3.to_excel(writer,"Sheet1",na_rep="",startrow=df_info.shape[0]+2)
df_ut_d2.to_excel(writer,"Sheet1",na_rep="",startrow=df_info.shape[0]+df_ut_d3.shape[0]+4)
df_vsp_d3.to_excel(writer,"Sheet1",na_rep="",startrow=df_info.shape[0]+df_ut_d3.shape[0]+df_ut_d2.shape[0]+6)
df_vsp_d2.to_excel(writer,"Sheet1",na_rep="",startrow=df_info.shape[0]+df_ut_d3.shape[0]+df_ut_d2.shape[0]+df_vsp_d3.shape[0]+8)
df_sample.to_excel(writer,"Sheet1",na_rep="",startrow=df_info.shape[0]+df_ut_d3.shape[0]+df_ut_d2.shape[0]+df_vsp_d3.shape[0]+df_vsp_d2.shape[0]+10)

df_cal.to_excel(writer,"Sheet1",na_rep="",startrow=3,startcol=df_ut_d3.shape[1]+2)
df_punch.to_excel(writer,"Sheet1",na_rep="",startrow=3+df_cal.shape[0]+2,startcol=df_ut_d3.shape[1]+2)
df_abvd.to_excel(writer,"Sheet1",na_rep="",startrow=3+df_cal.shape[0]+df_punch.shape[0]+4,startcol=df_ut_d3.shape[1]+2)

workbk=writer.book
worksh=writer.sheets['Sheet1']
chart1=workbk.add_chart({'type':'scatter'})
chart1.add_series({'name':['Sheet1',3,12],'categories':['Sheet1',4,11,11,11],'values':['Sheet1',4,12,11,12],'trendline':{'type':'linear','intercept':0,'display_equation':True,'display_r_squared':True}})
chart1.set_title({'name':'D3 Linearity'})
chart1.set_style(11)
worksh.insert_chart('L28',chart1)
chart2=workbk.add_chart({'type':'scatter'})
## Create chart2 with different notation
chart2.add_series({'name':'=Sheet1!$N$4','categories':'=Sheet1!$L$5:$L$12','values':'=Sheet1!$N$5:$N$12','trendline':{'type':'linear','intercept':0,'display_equation':True,'display_r_squared':True}})
chart2.set_title({'name':'D2 Linearity'})
chart2.set_style(11)
worksh.insert_chart('L43',chart2)

writer.save()

""" End of scripts """
