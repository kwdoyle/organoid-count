import os
import sys
import pandas as pd
import skimage.io as io
import colocalized_cell_count_functions as count

# TODO add RAGE to this, and the hopx with Rab (whatever that is)

def subset_img(img):
    y_step = img.shape[0] // 9 
    x_step = img.shape[1] // 4

    sections = []
    for i in range(9):
        for j in range(4):
            section = img[i*y_step:(i+1)*y_step, j*x_step:(j+1)*x_step]
            sections.append(section)

    return sections
  
  
def process_file(fl_lst):
    fl = fl_lst[0]
    img = io.imread(basedir + '/' + fl)
    img2 = (lambda x: x[:,:,:3] if x.shape[2] == 4 else x)(img)  # remove alpha channel if it exists
    # img_split = subset_img(img2)  # shouldn't need to subset these smaller images

    return img2  # img_split


def save_section_figs(section_img, section_name, section_num):
    if section_img is not None:
        count.savefig(section_img, savedir + "/sections/" + section_name + "/" + section_name + "_" + str(section_num) + ".png")
    else:
        print(f"No {section_name} image to save for section {section_num}")



basedir = os.path.normpath(sys.argv[1])
savedir = os.path.normpath(sys.argv[2])
toprocess = os.path.normpath(sys.argv[3])

min_clean_size = 5
close_radius = 3
min_watershed_dist = 15
rage_watershed = 7
dapi_watershed = 7
hopx_watershed = 8
aqp5_watershed = 8  #15 #8
dapi_hopx_watershed = 8
dapi_aqp5_watershed = 8
spc_watershed = 15
spc_cleansize = 1  # 5  # 1
sigma = 0.5
hopx_sigma = 0  # 0.5
spc_sigma = 0  # 0.2
aqp5_sigma = 0.5
spc_dilate = 0  # 1  # 10
font = 15
dapi_colocalize = False


# pull out the ID (the last folder name in the basedir argument) and use it in the output file name
id_nm = os.path.normpath(basedir).split(os.sep)[-1]
flnm = '/' + id_nm + '_' + 'count_output.xlsx'

# if os.path.exists(savedir + flnm):
#     print("File already analyzed")
#     sys.exit()
    
fls = os.listdir(basedir)


if toprocess == "rfp":
  os.makedirs(savedir + "/sections/", exist_ok=True)
  os.makedirs(savedir + "/out_sep/", exist_ok=True)
  os.makedirs(savedir + "/sections/erdr1/", exist_ok=True)
  os.makedirs(savedir + "/sections/tdt/", exist_ok=True)
  os.makedirs(savedir + "/sections/dapi/", exist_ok=True)
  os.makedirs(savedir + "/sections/rfp670/", exist_ok=True)
  os.makedirs(savedir + "/sections/full/", exist_ok=True)
  
  dapi_fl = [fl for fl in fls if "dapi" in fl.lower() and all(x.lower() not in fl.lower() for x in ['erdr1', 'tdt', 'spc', 'merging', 'merged'])]
  erdr1_fl = [fl for fl in fls if "erdr1" in fl.lower() and all(x.lower() not in fl.lower() for x in ['dapi', 'tdt', 'spc', 'merging', 'merged'])]
  tdt_fl = [fl for fl in fls if "rosa" in fl.lower() and all(x.lower() not in fl.lower() for x in ['dapi', 'erdr1', 'merging', 'merged'])]
  rfp670_fl = [fl for fl in fls if "rfp670" in fl.lower() and all(x.lower() not in fl.lower() for x in ['dapi', 'erdr1', 'tdt', 'merging', 'merged'])]
  full_fl = [fl for fl in fls if "merg" in fl.lower()]
  
  print("DAPI: " + str(dapi_fl))
  print("erdr1: " + str(erdr1_fl))
  print("tdt: " + str(tdt_fl))
  print("rfp670: " + str(rfp670_fl))
  print("full: " + str(full_fl))
  
  if len(dapi_fl) > 0:
      dapi_img_split = process_file(dapi_fl)
  else:
      dapi_img_split = []
  if len(erdr1_fl) > 0:
      erdr1_img_split = process_file(erdr1_fl)
  else:
      erdr1_img_split = []
  
  tdt_img_split = process_file(tdt_fl)
  rfp_img_split = process_file(rfp670_fl)
  full_img_split = process_file(full_fl)
  
  section_dat = pd.DataFrame()
  for i in range(0, len(rfp_img_split)):
      print("Processing section " + str(i))
      if len(erdr1_img_split) > 0:
          erdr1_sec = erdr1_img_split[i]
      else:
          erdr1_sec = None  # if no erdr1 image, set to None
      if len(dapi_img_split) > 0:
          dapi_sec = dapi_img_split[i]
      else:
          dapi_sec = None  # if no dapi image, set to None
      tdt_sec = tdt_img_split[i]
      rfp670_sec = rfp_img_split[i]
      full_sec = full_img_split[i]
  
      section = i+1
      try:
          erdr1_out = count.main(rfp670_img=rfp670_sec, tdt_img=tdt_sec, dapi_img=dapi_sec, erdr1_img=erdr1_sec, to_process="erdr1",
                              min_clean_size=min_clean_size, close_radius=close_radius, spc_watershed=spc_watershed, spc_sigma=spc_sigma, erdr1_sigma=erdr1_sigma,
                              min_watershed_dist=min_watershed_dist, dapi_watershed=dapi_watershed, spc_dilate=spc_dilate, spc_cleansize=spc_cleansize, sigma=sigma,
                              dapi_colocalize=False,
                              plot=True, savedir=savedir, section=i+1, font=font, stopcheck=False)
  
      except ValueError as e:
          print(e)
          continue
  
  
      alldat = erdr1_out
      alldat['section'] = i + 1
  
      section_dat = pd.concat([section_dat, alldat], ignore_index=True)
  
      # save sections for reference
      save_section_figs(erdr1_sec, "erdr1", section)
      save_section_figs(tdt_sec, "tdt", section)
      save_section_figs(dapi_sec, "dapi", section)
      save_section_figs(rfp670_sec, "rfp670", section)
      save_section_figs(full_sec, "full", section)
  
  
      flnmtmp = '/' + id_nm + '_' + 'sec_' + str(section) + '_count_output.xlsx'
      alldat.to_excel(savedir + '/out_sep/' + flnmtmp, index=False)
  
elif toprocess == "spc":
    dapi_fl = [fl for fl in fls if "dapi" in fl.lower() and all(x.lower() not in fl.lower() for x in ['spc', 'tdt', 'merging', 'merged'])]
    spc_fl = [fl for fl in fls if "spc" in fl.lower() and all(x.lower() not in fl.lower() for x in ['dapi', 'tdt', 'rosa', 'merging', 'merged'])]
    tdt_fl = [fl for fl in fls if "rosa" in fl.lower() and all(x.lower() not in fl.lower() for x in ['dapi', 'merging', 'merged'])]
    full_fl = [fl for fl in fls if "merg" in fl.lower()]
    
    print("DAPI: " + str(dapi_fl))
    print("spc: " + str(spc_fl))
    print("tdt: " + str(tdt_fl))
    print("full: " + str(full_fl))
    
    dapi_img = process_file(dapi_fl)
    spc_img = process_file(spc_fl)
    tdt_img = process_file(tdt_fl)
    full_img = process_file(full_fl)
    
    try:
        spc_out = count.main(spc_img=spc_img, tdt_img=tdt_img, dapi_img=dapi_img, to_process="spc",
                            min_clean_size=min_clean_size, close_radius=close_radius, spc_watershed=spc_watershed, spc_sigma=spc_sigma,
                            min_watershed_dist=min_watershed_dist, dapi_watershed=dapi_watershed, spc_dilate=spc_dilate, spc_cleansize=spc_cleansize, sigma=sigma,
                            dapi_colocalize=True,
                            plot=True, savedir=savedir, section=None, font=font)
    except ValueError as e:
        print(e)
    alldat = spc_out

    count.savefig(spc_img, savedir + "spc.png")
    count.savefig(tdt_img, savedir + "tdt.png")
    count.savefig(dapi_img, savedir + "dapi.png")
    count.savefig(full_img, savedir + "full.png")

elif toprocess == "lcn2":
  os.makedirs(savedir + "/sections/", exist_ok=True)
  os.makedirs(savedir + "/out_sep/", exist_ok=True)
  os.makedirs(savedir + "/sections/lcn2/", exist_ok=True)
  os.makedirs(savedir + "/sections/tdt/", exist_ok=True)
  os.makedirs(savedir + "/sections/dapi/", exist_ok=True)
  os.makedirs(savedir + "/sections/spc/", exist_ok=True)
  os.makedirs(savedir + "/sections/full/", exist_ok=True)

  dapi_fl = [fl for fl in fls if "dapi" in fl.lower() and all(x.lower() not in fl.lower() for x in ['hopx', 'tdt', 'spc', 'merging', 'merged'])]
  lcn2_fl = [fl for fl in fls if "lcn2" in fl.lower() and all(x.lower() not in fl.lower() for x in ['dapi', 'tdt', 'spc', 'merging', 'merged'])]
  tdt_fl = [fl for fl in fls if "rosa" in fl.lower() and all(x.lower() not in fl.lower() for x in ['dapi', 'hopx', 'merging', 'merged'])]
  spc_fl = [fl for fl in fls if "spc" in fl.lower() and all(x.lower() not in fl.lower() for x in ['dapi', 'hopx', 'tdt', 'rosa', 'merging', 'merged'])]
  full_fl = [fl for fl in fls if "merg" in fl.lower()]
  
  print("DAPI: " + str(dapi_fl))
  print("Lcn2: " + str(lcn2_fl))
  print("tdt: " + str(tdt_fl))
  print("spc: " + str(spc_fl))
  print("full: " + str(full_fl))
  
  lcn2_img_split = process_file(lcn2_fl)
  tdt_img_split = process_file(tdt_fl)
  dapi_img_split = process_file(dapi_fl)
  spc_img_split = process_file(spc_fl)
  full_img_split = process_file(full_fl)
  
  section_dat = pd.DataFrame()
  for i in range(0, len(lcn2_img_split)):
      print("Processing section " + str(i))
      lcn2_sec = lcn2_img_split[i]
      tdt_sec  = tdt_img_split[i]
      dapi_sec = dapi_img_split[i]
      spc_sec  = spc_img_split[i]
      full_sec = full_img_split[i]
  
      section = i+1
      try:
          lcn2_out = count.main(spc_img=spc_sec, tdt_img=tdt_sec, dapi_img=dapi_sec, lcn2_img=lcn2_sec, to_process="lcn2",
                              min_clean_size=min_clean_size, close_radius=close_radius, spc_watershed=spc_watershed, spc_sigma=spc_sigma,
                              min_watershed_dist=min_watershed_dist, dapi_watershed=dapi_watershed, spc_dilate=spc_dilate, spc_cleansize=spc_cleansize, sigma=sigma,
                              dapi_colocalize=True,
                              plot=True, savedir=savedir, section=i+1, font=font)

      except ValueError as e:
          print(e)
          continue
  
      alldat = lcn2_out
      alldat['section'] = i + 1
  
      section_dat = pd.concat([section_dat, alldat], ignore_index=True)
  
      # save sections for reference
      count.savefig(lcn2_sec, savedir + "/sections/lcn2/" + "lcn2_" + str(section) + ".png")
      count.savefig(tdt_sec, savedir + "/sections/tdt/" + "tdt_" + str(section) + ".png")
      count.savefig(dapi_sec, savedir + "/sections/dapi/" + "dapi_" + str(section) + ".png")
      count.savefig(spc_sec, savedir + "/sections/spc/" + "spc_" + str(section) + ".png")
      count.savefig(full_sec, savedir + "/sections/full/" + "full_" + str(section) + ".png")
  
      flnmtmp = '/' + id_nm + '_' + 'sec_' + str(section) + '_count_output.xlsx'
      alldat.to_excel(savedir + '/out_sep/' + flnmtmp, index=False)

elif toprocess == "hopx":
    dapi_fl = [fl for fl in fls if "dapi" in fl.lower() and all(x.lower() not in fl.lower() for x in ['spc', 'tdt', 'merging', 'merged'])]
    spc_fl = [fl for fl in fls if "spc" in fl.lower() and all(x.lower() not in fl.lower() for x in ['dapi', 'tdt', 'rosa', 'merging', 'merged'])]
    tdt_fl = [fl for fl in fls if "rosa" in fl.lower() and all(x.lower() not in fl.lower() for x in ['dapi', 'merging', 'merged'])]
    hopx_fl = [fl for fl in fls if "hopx" in fl.lower() and all(x.lower() not in fl.lower() for x in ['dapi', 'tdt', 'spc', 'merging', 'merged', 'rab'])]
    rab_fl = [fl for fl in fls if "rab" in fl.lower() and all(x.lower() not in fl.lower() for x in ['dapi', 'tdt', 'spc', 'merging', 'merged'])]
    full_fl = [fl for fl in fls if "merg" in fl.lower()]
    print("DAPI: " + str(dapi_fl))
    print("spc: " + str(spc_fl))
    print("tdt: " + str(tdt_fl))
    print("hopx: " + str(hopx_fl))
    print("rab: " + str(rab_fl))
    print("full: " + str(full_fl))
    
    dapi_img = process_file(dapi_fl)
    hopx_img = process_file(hopx_fl)
    full_img = process_file(full_fl)

    if len(spc_fl) > 0:
        spc_img = process_file(spc_fl)
    else:
        spc_img = None
    if len(tdt_fl) > 0:
        tdt_img = process_file(tdt_fl)
    else:
        tdt_img = None
    if len(rab_fl) > 0:
        rab_img = process_file(rab_fl)
    else:
        rab_img = None
    try:
        hopx_out = count.main(spc_img=spc_img, tdt_img=tdt_img, dapi_img=dapi_img, hopx_img=hopx_img, rab_img=rab_img, to_process="hopx",
                            min_clean_size=min_clean_size, close_radius=close_radius, spc_watershed=spc_watershed, spc_sigma=spc_sigma, 
                            hopx_watershed=hopx_watershed, dapi_hopx_watershed=dapi_hopx_watershed,
                            min_watershed_dist=min_watershed_dist, dapi_watershed=dapi_watershed, spc_dilate=spc_dilate, spc_cleansize=spc_cleansize, sigma=sigma, hopx_sigma=hopx_sigma,
                            dapi_colocalize=dapi_colocalize,
                            plot=True, savedir=savedir, section=None, font=font)
    except ValueError as e:
        print(e)

    alldat = hopx_out
    count.savefig(dapi_img, savedir + "dapi.png")
    count.savefig(hopx_img, savedir + "hopx.png")
    count.savefig(full_img, savedir + "full.png")
    if len(spc_fl) > 0:
        count.savefig(spc_img, savedir + "spc.png")
    if len(tdt_fl) > 0:
        count.savefig(tdt_img, savedir + "tdt.png")
    if len(rab_fl) > 0:
        count.savefig(rab_img, savedir + "rab.png")

elif toprocess == "aqp5":
    dapi_fl = [fl for fl in fls if "dapi" in fl.lower() and all(x.lower() not in fl.lower() for x in ['tdt', 'aqp5', 'merging', 'merged'])]
    aqp5_fl = [fl for fl in fls if "aqp5" in fl.lower() and all(x.lower() not in fl.lower() for x in ['dapi', 'tdt', 'merging', 'merged'])]
    tdt_fl = [fl for fl in fls if "rosa" in fl.lower() and all(x.lower() not in fl.lower() for x in ['dapi', 'aqp5', 'merging', 'merged'])]
    spc_fl = [fl for fl in fls if "spc" in fl.lower() and all(x.lower() not in fl.lower() for x in ['dapi', 'tdt', 'rosa', 'merging', 'merged'])]
    full_fl = [fl for fl in fls if "merg" in fl.lower()]
    
    print("DAPI: " + str(dapi_fl))
    print("Aqp5: " + str(aqp5_fl))
    print("tdt: " + str(tdt_fl))
    print("spc: " + str(spc_fl))
    print("full: " + str(full_fl))
    
    dapi_img = process_file(dapi_fl)
    aqp5_img = process_file(aqp5_fl)
    tdt_img = process_file(tdt_fl)
    full_img = process_file(full_fl)
    if len(spc_fl) > 0:
        spc_img = process_file(spc_fl)
    else:
        spc_img = None
    
    try:
        aqp5_out = count.main(aqp5_img=aqp5_img, tdt_img=tdt_img, spc_img=spc_img, dapi_img=dapi_img, to_process="aqp5",
                            min_clean_size=min_clean_size, close_radius=close_radius,
                            min_watershed_dist=min_watershed_dist, aqp5_watershed=aqp5_watershed, dapi_watershed=dapi_watershed, dapi_aqp5_watershed=dapi_aqp5_watershed,
                            sigma=sigma, aqp5_sigma=aqp5_sigma,
                            spc_watershed=spc_watershed, spc_sigma=spc_sigma, spc_cleansize=spc_cleansize, spc_dilate=spc_dilate,
                            dapi_colocalize=dapi_colocalize,
                            plot=True, savedir=savedir, section=None, font=font)

    except ValueError as e:
        print(e)
    alldat = aqp5_out
    
    count.savefig(aqp5_img, savedir + "aqp5.png")
    count.savefig(dapi_img, savedir + "dapi.png")
    count.savefig(full_img, savedir + "full.png")
    if len(spc_fl) > 0:
        count.savefig(spc_img, savedir + "spc.png")
    if len(tdt_fl) > 0:
        count.savefig(tdt_img, savedir + "tdt.png")

elif toprocess == "rage":
    dapi_fl = [fl for fl in fls if "dapi" in fl.lower() and all(x.lower() not in fl.lower() for x in ['rage', 'tdt', 'merging', 'merged'])]
    rage_fl = [fl for fl in fls if "rage" in fl.lower() and all(x.lower() not in fl.lower() for x in ['dapi', 'tdt', 'merging', 'merged'])]
    tdt_fl = [fl for fl in fls if "rosa" in fl.lower() and all(x.lower() not in fl.lower() for x in ['dapi', 'rage', 'merging', 'merged'])]
    spc_fl = [fl for fl in fls if "spc" in fl.lower() and all(x.lower() not in fl.lower() for x in ['dapi', 'tdt', 'rosa', 'rage', 'merging', 'merged'])]
    full_fl = [fl for fl in fls if "merg" in fl.lower()]
    
    print("DAPI: " + str(dapi_fl))
    print("Rage: " + str(rage_fl))
    print("tdt: " + str(tdt_fl))
    print("spc: " + str(spc_fl))
    print("full: " + str(full_fl))
    
    dapi_img = process_file(dapi_fl)
    rage_img = process_file(rage_fl)
    tdt_img = process_file(tdt_fl)
    full_img = process_file(full_fl)
    if len(spc_fl) > 0:
        spc_img = process_file(spc_fl)
    else:
        spc_img = None
    
    try:
        rage_out = count.main(rage_img=rage_img, tdt_img=tdt_img, spc_img=spc_img, dapi_img=dapi_img, to_process="rage",
                            min_clean_size=min_clean_size, close_radius=close_radius,
                            min_watershed_dist=min_watershed_dist, rage_watershed=rage_watershed, dapi_watershed=dapi_watershed,
                            sigma=sigma,
                            spc_watershed=spc_watershed, spc_sigma=spc_sigma, spc_cleansize=spc_cleansize, spc_dilate=spc_dilate,
                            dapi_colocalize=dapi_colocalize,
                            plot=True, savedir=savedir, section=None, font=font)
    except ValueError as e:
        print(e)
    alldat = rage_out

    count.savefig(rage_img, savedir + "rage.png")
    count.savefig(dapi_img, savedir + "dapi.png")
    count.savefig(full_img, savedir + "full.png")
    count.savefig(tdt_img, savedir + "tdt.png")
    if len(spc_fl) > 0:
        count.savefig(spc_img, savedir + "spc.png")

else:
  raise ValueError(f"Invalid 'toprocess' argument: {toprocess}")


## put image section number as first column and save table
# secnum = section_dat.pop('section')
# section_dat.insert(0, 'section', secnum)
# save section data
#section_dat.to_excel(savedir + flnm, index=False)

alldat.to_excel(savedir + flnm, index=False)

