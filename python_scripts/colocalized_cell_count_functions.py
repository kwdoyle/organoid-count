import os
os.environ["OMP_NUM_THREADS"] = '4'
import gc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from copy import deepcopy
from scipy import ndimage as ndi
from skimage.morphology import disk
from skimage.measure import regionprops, label
from skimage import color, feature, filters, measure, morphology, segmentation
from skimage.filters import threshold_multiotsu
from skimage import color


def plot_img(img):
    plt.imshow(img  , cmap='Greys_r')
    plt.show()


def savefig(img, fname):
    fig = plt.figure(frameon=False)

    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)

    ax.imshow(img, aspect='auto')

    fig.savefig(fname, dpi=1200)

    fig.clf()
    plt.close()


def make_colocalized_mask(img1, img2):
    # where the two masks overlap
    overlap = img1 * img2

    return overlap


def clean_mask(msk, minsize=64):
    tmp = morphology.remove_small_objects(measure.label(msk), min_size=minsize)
    tmp[np.nonzero(tmp)] = 1
    # convert back to int32
    tmp = tmp.astype(np.int32, copy=False)

    return tmp


def blur_mask(msk, sigma=1, return_all_1=True):
    msk_blur = filters.gaussian(msk, sigma=sigma)
    if return_all_1:
        msk_blur[np.nonzero(msk_blur)] = 1
        msk_blur = msk_blur.astype(np.int32, copy=False)

    return msk_blur


def make_watershed(msk, savedir, min_dist, plot=True, plt_title=None, show_numbers=True, interactive=False, font=0.005):
    if min_dist is None:
        raise ValueError("min_dist must be provided for watershed segmentation.")
    
    dist = ndi.distance_transform_edt(msk)
    coords = feature.peak_local_max(dist, min_distance=min_dist, footprint=np.ones((3, 3)), labels=msk)
    mask1 = np.zeros(dist.shape, dtype=bool)
    mask1[tuple(coords.T)] = True
    markers, _ = ndi.label(mask1)
    labels = segmentation.watershed(-dist, markers, mask=msk)
    labeled_cells = label(labels)

    if plot:
        fig, ax = plt.subplots(figsize=(9.60, 5.40))  # adjust the size as needed

        ax.imshow(msk, cmap=plt.cm.gray)  # display the mask
        ax.set_title(plt_title + ' n cells: ' + str(labels.max()))

        if show_numbers:
            # annotate each cell with its label number
            for region in regionprops(labeled_cells):
                y, x = region.centroid
                ax.text(x, y-15, str("."), color='red', fontsize=font, ha='center', va='center')

        ax.set_axis_off()

        fig.tight_layout()
        plt.savefig(savedir + '/' + plt_title + '.png', dpi=1200)
        if interactive:
            plt.show()
        fig.clf()
        plt.close('all')
        gc.collect()

    del labeled_cells

    return labels


def dapi(dapi_msk, plot, savedir, min_clean_size, min_watershed_dist, sigma, section=None, font=0.005, interactive=False):
    # just clean and count raw dapi cells
    dapi_msk_clean = blur_mask(clean_mask(dapi_msk, minsize=min_clean_size), sigma=sigma)
    if section is not None:
        title = "dapi_" + str(section)
    else:
        title = "dapi"
        
    n_dapi = make_watershed(msk=dapi_msk_clean, min_dist=min_watershed_dist, plot=plot, plt_title=title, savedir=savedir, interactive=interactive, font=font).max()

    return  n_dapi, dapi_msk_clean


def spc(spc_msk, plot, savedir, min_clean_size, dilate_radius, min_watershed_dist, sigma, clean=False, section=None, font=0.005, title="spc"):
    # always clean first
    spc_msk_clean = clean_mask(msk=spc_msk, minsize=min_clean_size)
    spc_msk_dilate = morphology.isotropic_dilation(spc_msk_clean, radius=dilate_radius)
    # test shrinking before and after bluring
    # bluring after seems better
    spc_msk_erode = morphology.binary_erosion(spc_msk_dilate, footprint=disk(3))
    spc_msk2 = blur_mask(msk=spc_msk_erode, sigma=sigma)

    n_spc = make_watershed(msk=spc_msk2, min_dist=min_watershed_dist, plot=plot, plt_title=title, savedir=savedir, font=font).max()

    return  n_spc, spc_msk2


def tdt(tdt_msk, plot, savedir, close_radius, min_clean_size, single_watershed, sigma, spc_watershed=None, spc_msk=None, 
        spc_sigma=None, spc_dilate=None, spc_cleansize=None, section=None, font=0.005, run_externally=True, run_spc=True):
  
    tdt_msk_clean = blur_mask(clean_mask(tdt_msk, minsize=min_clean_size), sigma=sigma)
    # tdt alone
    if section is not None:
        title = "tdt_" + str(section)
    else:
        title = "tdt"
    
    n_tdt = make_watershed(msk=tdt_msk_clean, min_dist=single_watershed, plot=plot, plt_title=title, savedir=savedir, font=font).max()
    # use function to get spc mask
    if run_spc:
        if section is not None:
            title = "spc_" + str(section)
        else:
            title = "spc"
        n_spc, spc_msk2 = spc(spc_msk=spc_msk, plot=plot, savedir=savedir, min_clean_size=spc_cleansize, 
                            min_watershed_dist=spc_watershed, dilate_radius=spc_dilate, sigma=spc_sigma,
                            section=section, font=font, title=title, clean=True)
     
        tdt_spc_msk = make_colocalized_mask(img1=tdt_msk_clean, img2=spc_msk2)
        tdt_spc_blur = blur_mask(msk=tdt_spc_msk, sigma=sigma)
      
        if section is not None:
            title = "tdt_spc_colocalized_" + str(section)
        else:
            title = "tdt_spc_colocalized"
            
        n_tdt_spc = make_watershed(msk=tdt_spc_blur, min_dist=single_watershed, plot=plot, plt_title=title, savedir=savedir, font=font, interactive=False).max()
        outdf = pd.DataFrame([{"n_tdt": n_tdt, "n_spc": n_spc, "n_tdt_spc": n_tdt_spc}])  #, "n_spc": n_spc}])

        if run_externally:
            return outdf
        else:
            return tdt_msk_clean, spc_msk2, n_tdt, n_spc, n_tdt_spc, tdt_spc_blur
    else:
        return tdt_msk_clean, n_tdt


def lcn2(lcn2_msk, tdt_msk, spc_msk, dapi_msk, plot, savedir, sigma, spc_sigma, close_radius, min_clean_size, single_watershed, dapi_watershed, spc_watershed, spc_cleansize, spc_dilate, section=None, font=0.005):
    print("processing dapi")
    n_dapi, dapi_msk_clean = dapi(dapi_msk=dapi_msk, plot=plot, savedir=savedir, sigma=sigma, font=font,
                                  min_watershed_dist=dapi_watershed, min_clean_size=min_clean_size, section=section)

    print("processing lcn2")
    lcn2_msk_clean = blur_mask(clean_mask(lcn2_msk, minsize=min_clean_size), sigma=sigma)

    if section is not None:
        title = "lcn2_" + str(section)
    else:
        title = "lcn2"
        
    n_lcn2 = make_watershed(msk=lcn2_msk_clean, min_dist=single_watershed, plot=plot, plt_title=title, savedir=savedir, font=font).max()

    # tdt
    print("processing tdt and spc")
    tdt_msk_clean, spc_msk2, n_tdt, n_spc, n_tdt_spc, tdt_spc_msk = tdt(tdt_msk=tdt_msk, spc_msk=spc_msk, plot=plot, savedir=savedir,
                                                            close_radius=close_radius, min_clean_size=min_clean_size, 
                                                            single_watershed=single_watershed, spc_watershed=spc_watershed,
                                                            sigma=sigma, spc_sigma=spc_sigma, spc_dilate=spc_dilate, 
                                                            spc_cleansize=spc_cleansize, section=section, font=font, run_externally=False)


    print("Processing lcn2 & tdt")
    lcn2_tdt_msk = blur_mask(make_colocalized_mask(img1=lcn2_msk_clean, img2=tdt_msk_clean), sigma=sigma)

    if section is not None:
        title = "lcn2_tdt_colocalized_" + str(section)
    else:
        title = "lcn2_tdt_colocalized_"
    n_lcn2_tdt = make_watershed(msk=lcn2_tdt_msk, min_dist=single_watershed, plot=plot, plt_title=title, savedir=savedir, font=font, interactive=False).max()


    print("Processing lcn2 & spc")
    lcn2_spc_msk = blur_mask(make_colocalized_mask(img1=lcn2_msk_clean, img2=spc_msk2), sigma=sigma)

    if section is not None:
        title = "lcn2_spc_colocalized_" + str(section)
    else:
        title = "lcn2_spc_colocalized_"
    n_lcn2_spc = make_watershed(msk=lcn2_spc_msk, min_dist=single_watershed, plot=plot, plt_title=title, savedir=savedir, font=font, interactive=False).max()


    print("Processing lcn2 & tdt & spc")
    lcn2_tdt_spc_msk = blur_mask(make_colocalized_mask(img1=lcn2_tdt_msk, img2=spc_msk2), sigma=sigma)

    if section is not None:
        title = "lcn2_tdt_spc_colocalized_" + str(section)
    else:
        title = "lcn2_tdt_spc_colocalized_"
    n_lcn2_tdt_spc = make_watershed(msk=lcn2_tdt_spc_msk, min_dist=single_watershed, plot=plot, plt_title=title, savedir=savedir, font=font, interactive=False).max()


    print("Processing (spc & tdt) - lcn2")
    # dilate lcn2 mask first
    lcn2_msk_dilate = morphology.isotropic_dilation(lcn2_msk_clean, radius=12)  # 3
    tdt_spc_minus_lcn2 = np.subtract(np.array(tdt_spc_msk, dtype="int16"), np.array(lcn2_msk_dilate, dtype="int16"))
    tdt_spc_minus_lcn2[np.where(tdt_spc_minus_lcn2 < 0)] = 0 
    tdt_spc_minus_lcn2_clean = clean_mask(tdt_spc_minus_lcn2, minsize=min_clean_size)

    # close and clean
    # still close after subtracting
    tdt_spc_minus_lcn2_closed = blur_mask(morphology.isotropic_closing(tdt_spc_minus_lcn2_clean, radius=close_radius), sigma=sigma)
    # dilate and erode
    tdt_spc_minus_lcn2_closed_dilate = morphology.isotropic_dilation(tdt_spc_minus_lcn2_closed, radius=12)
    tdt_spc_minus_lcn2_closed_erode = morphology.binary_erosion(tdt_spc_minus_lcn2_closed_dilate, footprint=disk(3))

    if section is not None:
        title = "tdt_spc_colocalized_minus_lcn2" + str(section)
    else:
        title = "tdt_spc_colocalized_minus_lcn2"
    n_tdt_spc_minus_lcn2 = make_watershed(msk=tdt_spc_minus_lcn2_closed_erode, min_dist=single_watershed, plot=plot, plt_title=title, savedir=savedir, font=font, interactive=False).max()

    outdf = pd.DataFrame([{"n_dapi": n_dapi, "n_lcn2": n_lcn2, "n_tdt": n_tdt, "n_spc": n_spc, "n_tdt_spc": n_tdt_spc, "n_lcn2_tdt": n_lcn2_tdt, "n_lcn2_spc": n_lcn2_spc,
                           "n_lcn2_tdt_spc": n_lcn2_tdt_spc, "n_tdt_spc_minus_lcn2": n_tdt_spc_minus_lcn2}])

    return outdf


def erdr1(erdr1_msk, tdt_msk, rfp670_msk, dapi_msk, plot, savedir, sigma, spc_sigma, erdr1_sigma, close_radius, min_clean_size, single_watershed, dapi_watershed, 
          spc_watershed, spc_cleansize, spc_dilate, spc_msk=None, dapi_colocalize=False, section=None, font=0.005):

    if dapi_msk:
        print("processing dapi")
        n_dapi, dapi_msk_clean = dapi(dapi_msk=dapi_msk, plot=plot, savedir=savedir, sigma=sigma, font=font,
                                    min_watershed_dist=dapi_watershed, min_clean_size=min_clean_size, section=section)
    else:
        n_dapi = 0

    if erdr1_msk:
        print("processing erdr1")
        erdr1_msk_clean = blur_mask(clean_mask(erdr1_msk, minsize=min_clean_size), sigma=erdr1_sigma)

        if section is not None:
            title = "erdr1_" + str(section)
        else:
            title = "erdr1"
            
        n_erdr1 = make_watershed(msk=erdr1_msk_clean, min_dist=single_watershed, plot=plot, plt_title=title, savedir=savedir, font=font).max()
    else:
        n_erdr1 = 0
    

    print("processing rfp670")
    rfp670_msk_clean = blur_mask(clean_mask(rfp670_msk, minsize=min_clean_size), sigma=erdr1_sigma)

    if section is not None:
        title = "rfp670_" + str(section)
    else:
        title = "rfp670"
        
    n_rfp670 = make_watershed(msk=rfp670_msk_clean, min_dist=single_watershed, plot=plot, plt_title=title, savedir=savedir, font=font).max()

    # tdt
    print("processing tdt")
    tdt_msk_clean, n_tdt = tdt(tdt_msk=tdt_msk, spc_msk=spc_msk, plot=plot, savedir=savedir,
                                                            close_radius=close_radius, min_clean_size=min_clean_size, 
                                                            single_watershed=single_watershed, spc_watershed=spc_watershed,
                                                            sigma=sigma, spc_sigma=spc_sigma, spc_dilate=spc_dilate, run_spc=False, 
                                                            spc_cleansize=spc_cleansize, section=section, font=font, run_externally=False)


    if dapi_msk and erdr1_msk:
        print("Processing erdr1 & dapi")
        erdr1_dapi_msk = make_colocalized_mask(img1=erdr1_msk_clean, img2=dapi_msk_clean)

        erdr1_dapi_msk_blur = blur_mask(erdr1_dapi_msk, sigma=erdr1_sigma)

        if section is not None:
            title = "erdr1_dapi_colocalized_" + str(section)
        else:
            title = "erdr1_dapi_colocalized"
            
        n_erdr1_dapi = make_watershed(msk=erdr1_dapi_msk_blur, min_dist=dapi_watershed, plot=plot, plt_title=title, savedir=savedir, font=font, interactive=False).max()
    else:
        n_erdr1_dapi = 0

    if dapi_colocalize:
        n_erdr1 = n_erdr1_dapi


    if erdr1_msk:
        print("Processing erdr1 & tdt")
        if dapi_colocalize:
            erdr1_tdt_msk = blur_mask(make_colocalized_mask(img1=erdr1_dapi_msk_blur, img2=tdt_msk_clean), sigma=sigma)
        else:
            erdr1_tdt_msk = blur_mask(make_colocalized_mask(img1=erdr1_msk_clean, img2=tdt_msk_clean), sigma=sigma)

        if section is not None:
            title = "erdr1_tdt_colocalized_" + str(section)
        else:
            title = "erdr1_tdt_colocalized"
            
        n_erdr1_tdt = make_watershed(msk=erdr1_tdt_msk, min_dist=single_watershed, plot=plot, plt_title=title, savedir=savedir, font=font, interactive=False).max()
    else:
        n_erdr1_tdt = 0


    if erdr1_msk:
        print("Processing erdr1 & RFP670")
        if dapi_colocalize:
            erdr1_rfp670_msk = blur_mask(make_colocalized_mask(img1=erdr1_dapi_msk_blur, img2=rfp670_msk_clean), sigma=sigma)
        else:
            erdr1_rfp670_msk = blur_mask(make_colocalized_mask(img1=erdr1_msk_clean, img2=rfp670_msk_clean), sigma=sigma)

        if section is not None:
            title = "erdr1_rfp670_colocalized_" + str(section)
        else:
            title = "erdr1_rfp670_colocalized"
            
        n_erdr1_rfp670 = make_watershed(msk=erdr1_rfp670_msk, min_dist=single_watershed, plot=plot, plt_title=title, savedir=savedir, font=font, interactive=False).max()    
    else:
        n_erdr1_rfp670 = 0
    

    print("Processing tdt & RFP670")
    tdt_rfp670_msk = blur_mask(make_colocalized_mask(img1=tdt_msk_clean, img2=rfp670_msk_clean), sigma=sigma)

    if section is not None:
        title = "tdt_rfp670_colocalized_" + str(section)
    else:
        title = "tdt_rfp670_colocalized"
        
    n_tdt_rfp670 = make_watershed(msk=tdt_rfp670_msk, min_dist=single_watershed, plot=plot, plt_title=title, savedir=savedir, font=font, interactive=False).max()    


    if erdr1_msk:
        print("Processing erdr1 & rfp670 & tdt")
        if dapi_colocalize:
            erdr1_tdt_rfp670_msk = blur_mask(make_colocalized_mask(img1=erdr1_dapi_msk_blur, img2=tdt_rfp670_msk), sigma=sigma)
        else:
            erdr1_tdt_rfp670_msk = blur_mask(make_colocalized_mask(img1=erdr1_msk_clean, img2=tdt_rfp670_msk), sigma=sigma)

        if section is not None:
            title = "erdr1_tdt_rfp670_colocalized_" + str(section)
        else:
            title = "erdr1_tdt_rfp670_colocalized"
            
        n_erdr1_tdt_rfp670 = make_watershed(msk=erdr1_tdt_rfp670_msk, min_dist=single_watershed, plot=plot, plt_title=title, savedir=savedir, font=font, interactive=False).max()
    else:
        n_erdr1_tdt_rfp670 = 0


    outdf = pd.DataFrame([{"n_dapi": n_dapi, "n_erdr1": n_erdr1, "n_tdt": n_tdt, "n_rfp670": n_rfp670, #"n_erdr1_dapi": n_erdr1_dapi, 
                           "n_erdr1_tdt": n_erdr1_tdt, "n_erdr1_rfp670": n_erdr1_rfp670, "n_tdt_rfp670": n_tdt_rfp670, "n_erdr1_tdt_rfp670": n_erdr1_tdt_rfp670}])

    return outdf


def hopx(hopx_msk, dapi_msk, tdt_msk, spc_msk, rab_msk, plot, savedir, sigma, close_radius, min_clean_size, spc_cleansize, single_watershed, dapi_watershed, dapi_hopx_watershed, 
         spc_watershed, spc_dilate, spc_sigma, dapi_sigma, hopx_spc_clean, hopx_tdt_clean, dapi_colocalize=False, font=0.005, section=None):  # spc_close=3
    ## this should give counts of raw hopx cells, hopx on dapi cells, hopix/dapi on tdt, hopix/dapi on tdt minus spc
    ## and give intersection coefs for each pair

    # first create cleaned dapi mask
    print("Processing dapi")
    n_dapi, dapi_msk_clean = dapi(dapi_msk=dapi_msk, plot=plot, savedir=savedir, sigma=dapi_sigma, min_watershed_dist=dapi_watershed,
                                  min_clean_size=min_clean_size, section=section, font=font)

    ## hopx alone:
    # clean hopx mask first
    print("processing hopx")
    hopx_closed = morphology.isotropic_closing(hopx_msk, radius=close_radius)
    hopx_msk_clean = blur_mask(clean_mask(hopx_closed, minsize=min_clean_size), sigma=sigma)
    # clean and check hopx independently as well
    if section is not None:
        title = "hopx_" + str(section)
    else:
        title = "hopx"
    n_hopx = make_watershed(msk=hopx_msk_clean, min_dist=single_watershed, plot=plot, plt_title=title, savedir=savedir, font=font).max()


    print("Processing dapi & hopx")
    hopx_dapi_msk = blur_mask(make_colocalized_mask(img1=dapi_msk_clean, img2=hopx_msk_clean), sigma=sigma)

    if section is not None:
        title = "hopx_dapi_colocalized_" + str(section)
    else:
        title = "hopx_dapi_colocalized"

    n_hopx_dapi = make_watershed(msk=hopx_dapi_msk, min_dist=dapi_hopx_watershed, plot=plot, plt_title=title, savedir=savedir, font=font, interactive=False).max()


    if tdt_msk is not None and spc_msk is not None:
        print("Processing tdt & spc")
        tdt_msk_clean, spc_msk2, n_tdt, n_spc, n_tdt_spc, tdt_spc_msk = tdt(tdt_msk=tdt_msk, spc_msk=spc_msk, plot=plot, savedir=savedir,
                                                                close_radius=close_radius, min_clean_size=min_clean_size, 
                                                                single_watershed=single_watershed, spc_watershed=spc_watershed,
                                                                sigma=sigma, spc_sigma=spc_sigma, spc_dilate=spc_dilate, 
                                                                spc_cleansize=spc_cleansize, section=section, font=font, run_externally=False)
        # hopx/dapi and tdt overlap
        if dapi_colocalize:
            hopx_tdt_msk = make_colocalized_mask(img1=hopx_dapi_msk, img2=tdt_msk_clean)
        else:
            hopx_tdt_msk = make_colocalized_mask(img1=hopx_msk_clean, img2=tdt_msk_clean)
        # if I need to close after subtracting spc, I will need to close here as well
        # in case some almost-connected spot here remains the same after the subtraction,
        # the closing after subtracting will now count it as one when it should still be two.
        # note: the watershed distance might make this moot anyway
        hopx_tdt_closed = blur_mask(morphology.isotropic_closing(hopx_tdt_msk, radius=close_radius), sigma=sigma)
        if section is not None:
            title = "hopx_tdt_colocalized_" + str(section)
        else:
            title = "hopx_tdt_colocalized"

        n_hopx_tdt = make_watershed(msk=hopx_tdt_closed, min_dist=single_watershed, plot=plot, plt_title=title, savedir=savedir, font=font, interactive=False).max()

        # hopx and tdt minus spc:
        print("Processing tdt - spc")
        hopx_tdt_minus_spc = np.subtract(np.array(hopx_tdt_closed, dtype="int16"), np.array(spc_msk2, dtype="int16"))
        # del hopx_tdt_closed_clean
        hopx_tdt_minus_spc[np.where(hopx_tdt_minus_spc < 0)] = 0 
        hopx_tdt_minus_spc = clean_mask(hopx_tdt_minus_spc, minsize=hopx_tdt_clean)

        # some other subtractions
        if dapi_colocalize:
            hopx_dapi_minus_spc = np.subtract(np.array(hopx_dapi_msk, dtype="int16"), np.array(spc_msk2, dtype="int16"))
        else:
            hopx_dapi_minus_spc = np.subtract(np.array(hopx_msk_clean, dtype="int16"), np.array(spc_msk2, dtype="int16"))
        hopx_dapi_minus_spc[np.where(hopx_dapi_minus_spc < 0)] = 0 
        hopx_dapi_minus_spc = clean_mask(hopx_dapi_minus_spc, minsize=hopx_spc_clean)

        # close and clean
        # still close after subtracting
        hopx_tdt_minus_spc_closed = blur_mask(morphology.isotropic_closing(hopx_tdt_minus_spc, radius=close_radius), sigma=sigma)
        hopx_dapi_minus_spc_closed = blur_mask(morphology.isotropic_closing(hopx_dapi_minus_spc, radius=close_radius), sigma=sigma)
        # hopx_tdt_minus_spc_closed_clean = clean_mask(hopx_tdt_minus_spc_closed, minsize=min_clean_size)
        # del hopx_tdt_minus_spc_closed, hopx_tdt_minus_spc
        if section is not None:
            title = "hopx_tdt_colocalized_minus_spc_" + str(section)
        else:
            title = "hopx_tdt_colocalized_minus_spc"
        # hopx_tdt_minus_spc_closed_watershed_labels = make_watershed(msk=hopx_tdt_minus_spc_closed, min_dist=min_watershed_dist, plot=plot, plt_title=title, savedir=savedir, interactive=False)
        # n_hopx_tdt_minus_spc = hopx_tdt_minus_spc_closed_watershed_labels.max()
        n_hopx_tdt_minus_spc = make_watershed(msk=hopx_tdt_minus_spc_closed, min_dist=single_watershed, plot=plot, plt_title=title, savedir=savedir, font=font, interactive=False).max()
        # n_hopx_tdt_minus_spc = hopx_tdt_minus_spc_closed_watershed_labels.max()
        # del hopx_tdt_minus_spc_closed_watershed_labels

        if section is not None:
            title = "hopx_dapi_colocalized_minus_spc_" + str(section)
        else:
            title = "hopx_dapi_colocalized_minus_spc"
        # hopx_tdt_minus_spc_closed_watershed_labels = make_watershed(msk=hopx_tdt_minus_spc_closed, min_dist=min_watershed_dist, plot=plot, plt_title=title, savedir=savedir, interactive=False)
        # n_hopx_tdt_minus_spc = hopx_tdt_minus_spc_closed_watershed_labels.max()
        n_hopx_dapi_minus_spc = make_watershed(msk=hopx_dapi_minus_spc_closed, min_dist=single_watershed, plot=plot, plt_title=title, savedir=savedir, font=font, interactive=False).max()
        

        # hopx dapi spc
        if dapi_colocalize:
            hopx_spc_msk = make_colocalized_mask(img1=hopx_dapi_msk, img2=spc_msk2)
        else:
            hopx_spc_msk = make_colocalized_mask(img1=hopx_msk_clean, img2=spc_msk2)
        # if I need to close after subtracting spc, I will need to close here as well
        # in case some almost-connected spot here remains the same after the subtraction,
        # the closing after subtracting will now count it as one when it should still be two.
        # note: the watershed distance might make this moot anyway
        hopx_spc_closed = blur_mask(morphology.isotropic_closing(hopx_spc_msk, radius=close_radius), sigma=sigma)
        if section is not None:
            title = "hopx_spc_colocalized_" + str(section)
        else:
            title = "hopx_spc_colocalized"
        # shed = make_watershed(msk=hopx_tdt_closed, min_dist=min_watershed_dist, plot=plot, plt_title=title, savedir=savedir, interactive=False)
        # n_hopx_tdt = shed.max()
        n_hopx_spc = make_watershed(msk=hopx_spc_closed, min_dist=single_watershed, plot=plot, plt_title=title, savedir=savedir, font=font, interactive=False).max()

        # hopx dapi tdt spc
        # use initial masks before closing and cleaning them--then close and clean this resulting mask
        hopx_tdt_spc_msk = make_colocalized_mask(img1=hopx_tdt_msk, img2=spc_msk2)
        # if I need to close after subtracting spc, I will need to close here as well
        # in case some almost-connected spot here remains the same after the subtraction,
        # the closing after subtracting will now count it as one when it should still be two.
        # note: the watershed distance might make this moot anyway
        hopx_tdt_spc_closed = blur_mask(morphology.isotropic_closing(hopx_tdt_spc_msk, radius=close_radius), sigma=sigma)
        if section is not None:
            title = "hopx_tdt_spc_colocalized_" + str(section)
        else:
            title = "hopx_tdt_spc_colocalized"
        # shed = make_watershed(msk=hopx_tdt_closed, min_dist=min_watershed_dist, plot=plot, plt_title=title, savedir=savedir, interactive=False)
        # n_hopx_tdt = shed.max()
        n_hopx_tdt_spc = make_watershed(msk=hopx_tdt_spc_closed, min_dist=single_watershed, plot=plot, plt_title=title, savedir=savedir, font=font, interactive=False).max()

    else:
        n_tdt = 0
        n_spc = 0
        n_tdt_spc = 0
        n_hopx_tdt = 0
        n_hopx_tdt_minus_spc = 0
        n_hopx_dapi_minus_spc = 0
        n_hopx_spc = 0
        n_hopx_tdt_spc = 0

    if rab_msk is not None:
        print("Processing hopx & rab")
        hopx_rab_msk = make_colocalized_mask(img1=hopx_msk_clean, img2=rab_msk) # TODO should use hopx + dapi mask?
        hopx_rab_closed = blur_mask(morphology.isotropic_closing(hopx_rab_msk, radius=close_radius), sigma=sigma)
        if section is not None:
            title = "hopx_rab_colocalized_" + str(section)
        else:
            title = "hopx_rab_colocalized"
        n_hopx_rab = make_watershed(msk=hopx_rab_closed, min_dist=single_watershed, plot=plot, plt_title=title, savedir=savedir, font=font, interactive=False).max()
    else:
        n_hopx_rab = 0

    # make df
    outdf = pd.DataFrame([{"n_dapi": n_dapi, "n_hopx": n_hopx, "n_spc": n_spc, "n_tdt": n_tdt, "n_tdt_spc": n_tdt_spc, "n_hopx_dapi": n_hopx_dapi, "n_hopx_dapi_tdt": n_hopx_tdt,
                           "n_hopx_dapi_tdt_spc": n_hopx_tdt_spc, "n_hopx_dapi_spc": n_hopx_spc, "n_hopx_dapi_tdt_minus_spc": n_hopx_tdt_minus_spc,
                           "n_hopx_dapi_minus_spc": n_hopx_dapi_minus_spc, "n_hopx_rab": n_hopx_rab}])

    return outdf


def aqp5(aqp5_msk, dapi_msk, tdt_msk, spc_msk, plot, savedir, sigma, aqp5_sigma, single_watershed, dapi_watershed, dapi_aqp5_watershed, aqp5_spc_clean, aqp5_tdt_clean, 
         spc_watershed, spc_sigma, spc_cleansize, spc_dilate, close_radius, min_clean_size, dapi_colocalize=False, section=None, font=0.005):
    ## new way, not counting per cell, but finding overall percentage of aqp5 area over other area

    # first create cleaned dapi mask
    print("processing dapi")
    n_dapi, dapi_msk_clean = dapi(dapi_msk=dapi_msk, plot=plot, savedir=savedir, sigma=sigma, font=font,
                                  min_watershed_dist=dapi_watershed, min_clean_size=min_clean_size, section=section)

    print("processing aqp5")
    aqp5_msk_clean = blur_mask(clean_mask(aqp5_msk, minsize=min_clean_size), sigma=aqp5_sigma)

    if section is not None:
        title = "aqp5_" + str(section)
    else:
        title = "aqp5"
    n_aqp5 = make_watershed(msk=aqp5_msk_clean, min_dist=single_watershed, plot=plot, plt_title=title, savedir=savedir, font=font).max()

    print("Processing aqp5 & dapi")
    aqp5_dapi_msk = blur_mask(make_colocalized_mask(img1=aqp5_msk_clean, img2=dapi_msk_clean), sigma=sigma)
    if section is not None:
        title = "aqp5_dapi_colocalized_" + str(section)
    else:
        title = "aqp5_dapi_colocalized"

    n_aqp5_dapi = make_watershed(msk=aqp5_dapi_msk, min_dist=dapi_aqp5_watershed, plot=plot, plt_title=title, savedir=savedir, font=font, interactive=False).max()
    
    # if have tdt or spc..
    if spc_msk is not None or tdt_msk is not None:
        # process both tdt and spc if both are given
        if spc_msk is not None and tdt_msk is not None:
            print("Processing tdt & spc")
            tdt_msk_clean, spc_msk2, n_tdt, n_spc, n_tdt_spc, tdt_spc_msk = tdt(tdt_msk=tdt_msk, spc_msk=spc_msk, plot=plot, savedir=savedir,
                                                                    close_radius=close_radius, min_clean_size=min_clean_size, 
                                                                    single_watershed=single_watershed, spc_watershed=spc_watershed,
                                                                    sigma=sigma, spc_sigma=spc_sigma, spc_dilate=spc_dilate, 
                                                                    spc_cleansize=spc_cleansize, section=section, font=font, run_externally=False)
        # else process only tdt if spc not given
        elif tdt_msk is not None and spc_msk is None:
            print("Processing tdt without spc")
            tdt_msk_clean, n_tdt = tdt(tdt_msk=tdt_msk, run_spc=False, plot=plot, savedir=savedir,
                                                                    close_radius=close_radius, min_clean_size=min_clean_size, 
                                                                    single_watershed=single_watershed, 
                                                                    sigma=sigma, 
                                                                    section=section, font=font, run_externally=False)

        # TODO put parts of these corresponding to whether spc is given or not
        # hopx/dapi and tdt overlap
        if dapi_colocalize:
            aqp5_tdt_msk = make_colocalized_mask(img1=aqp5_dapi_msk, img2=tdt_msk_clean)
        else:
            aqp5_tdt_msk = make_colocalized_mask(img1=aqp5_msk_clean, img2=tdt_msk_clean)
        # if I need to close after subtracting spc, I will need to close here as well
        # in case some almost-connected spot here remains the same after the subtraction,
        # the closing after subtracting will now count it as one when it should still be two.
        # note: the watershed distance might make this moot anyway
        aqp5_tdt_closed = blur_mask(morphology.isotropic_closing(aqp5_tdt_msk, radius=close_radius), sigma=sigma)
        if section is not None:
            title = "aqp5_tdt_colocalized_" + str(section)
        else:
            title = "aqp5_tdt_colocalized"

        n_aqp5_tdt = make_watershed(msk=aqp5_tdt_closed, min_dist=single_watershed, plot=plot, plt_title=title, savedir=savedir, font=font, interactive=False).max()

        if spc_msk is not None:
            # hopx and tdt minus spc:
            print("Processing tdt - spc")
            aqp5_tdt_minus_spc = np.subtract(np.array(aqp5_tdt_closed, dtype="int16"), np.array(spc_msk2, dtype="int16"))
            # del hopx_tdt_closed_clean
            aqp5_tdt_minus_spc[np.where(aqp5_tdt_minus_spc < 0)] = 0 
            aqp5_tdt_minus_spc = clean_mask(aqp5_tdt_minus_spc, minsize=aqp5_tdt_clean)

            aqp5_tdt_minus_spc_closed = blur_mask(morphology.isotropic_closing(aqp5_tdt_minus_spc, radius=close_radius), sigma=sigma)
            if section is not None:
                title = "aqp5_tdt_colocalized_minus_spc_" + str(section)
            else:
                title = "aqp5_tdt_colocalized_minus_spc"
            n_aqp5_tdt_minus_spc = make_watershed(msk=aqp5_tdt_minus_spc_closed, min_dist=single_watershed, plot=plot, plt_title=title, savedir=savedir, font=font, interactive=False).max()

            # some other subtractions
            if dapi_colocalize:
                aqp5_dapi_minus_spc = np.subtract(np.array(aqp5_dapi_msk, dtype="int16"), np.array(spc_msk2, dtype="int16"))
            else:
                aqp5_dapi_minus_spc = np.subtract(np.array(aqp5_msk_clean, dtype="int16"), np.array(spc_msk2, dtype="int16"))
            aqp5_dapi_minus_spc[np.where(aqp5_dapi_minus_spc < 0)] = 0 
            aqp5_dapi_minus_spc = clean_mask(aqp5_dapi_minus_spc, minsize=aqp5_spc_clean)

            # close and clean
            # still close after subtracting
            aqp5_dapi_minus_spc_closed = blur_mask(morphology.isotropic_closing(aqp5_dapi_minus_spc, radius=close_radius), sigma=sigma)
            # hopx_tdt_minus_spc_closed_clean = clean_mask(hopx_tdt_minus_spc_closed, minsize=min_clean_size)
            # del hopx_tdt_minus_spc_closed, hopx_tdt_minus_spc

            # n_hopx_tdt_minus_spc = hopx_tdt_minus_spc_closed_watershed_labels.max()
            # del hopx_tdt_minus_spc_closed_watershed_labels

            if section is not None:
                title = "aqp5_dapi_colocalized_minus_spc_" + str(section)
            else:
                title = "aqp5_dapi_colocalized_minus_spc"
            # hopx_tdt_minus_spc_closed_watershed_labels = make_watershed(msk=hopx_tdt_minus_spc_closed, min_dist=min_watershed_dist, plot=plot, plt_title=title, savedir=savedir, interactive=False)
            # n_hopx_tdt_minus_spc = hopx_tdt_minus_spc_closed_watershed_labels.max()
            n_aqp5_dapi_minus_spc = make_watershed(msk=aqp5_dapi_minus_spc_closed, min_dist=single_watershed, plot=plot, plt_title=title, savedir=savedir, font=font, interactive=False).max()
        

            # hopx dapi spc
            if dapi_colocalize:
                aqp5_spc_msk = make_colocalized_mask(img1=aqp5_dapi_msk, img2=spc_msk2)
            else:
                aqp5_spc_msk = make_colocalized_mask(img1=aqp5_msk_clean, img2=spc_msk2)
            # if I need to close after subtracting spc, I will need to close here as well
            # in case some almost-connected spot here remains the same after the subtraction,
            # the closing after subtracting will now count it as one when it should still be two.
            # note: the watershed distance might make this moot anyway
            aqp5_spc_closed = blur_mask(morphology.isotropic_closing(aqp5_spc_msk, radius=close_radius), sigma=sigma)
            if section is not None:
                title = "aqp5_spc_colocalized_" + str(section)
            else:
                title = "aqp5_spc_colocalized"
            # shed = make_watershed(msk=hopx_tdt_closed, min_dist=min_watershed_dist, plot=plot, plt_title=title, savedir=savedir, interactive=False)
            # n_hopx_tdt = shed.max()
            n_aqp5_spc = make_watershed(msk=aqp5_spc_closed, min_dist=single_watershed, plot=plot, plt_title=title, savedir=savedir, font=font, interactive=False).max()

            # hopx dapi tdt spc
            # use initial masks before closing and cleaning them--then close and clean this resulting mask
            aqp5_tdt_spc_msk = make_colocalized_mask(img1=aqp5_tdt_msk, img2=spc_msk2)
            # if I need to close after subtracting spc, I will need to close here as well
            # in case some almost-connected spot here remains the same after the subtraction,
            # the closing after subtracting will now count it as one when it should still be two.
            # note: the watershed distance might make this moot anyway
            aqp5_tdt_spc_closed = blur_mask(morphology.isotropic_closing(aqp5_tdt_spc_msk, radius=close_radius), sigma=sigma)
            if section is not None:
                title = "aqp5_tdt_spc_colocalized_" + str(section)
            else:
                title = "aqp5_tdt_spc_colocalized"
            n_aqp5_tdt_spc = make_watershed(msk=aqp5_tdt_spc_closed, min_dist=single_watershed, plot=plot, plt_title=title, savedir=savedir, font=font, interactive=False).max()

        else:
            n_spc = 0
            n_tdt_spc = 0
            n_aqp5_tdt_minus_spc = 0
            n_aqp5_dapi_minus_spc = 0
            n_aqp5_spc = 0
            n_aqp5_tdt_spc = 0

    else:
        n_tdt = 0
        n_spc = 0
        n_tdt_spc = 0
        n_aqp5_tdt = 0
        n_aqp5_tdt_minus_spc = 0
        n_aqp5_dapi_minus_spc = 0
        n_aqp5_spc = 0
        n_aqp5_tdt_spc = 0


    # make df
    outdf = pd.DataFrame([{"n_dapi": n_dapi, "n_aqp5": n_aqp5, "n_spc": n_spc, "n_tdt": n_tdt, "n_tdt_spc": n_tdt_spc, "n_aqp5_dapi": n_aqp5_dapi, "n_aqp5_dapi_tdt": n_aqp5_tdt,
                           "n_aqp5_dapi_tdt_spc": n_aqp5_tdt_spc, "n_aqp5_dapi_spc": n_aqp5_spc, "n_aqp5_dapi_tdt_minus_spc": n_aqp5_tdt_minus_spc,
                           "n_aqp5_dapi_minus_spc": n_aqp5_dapi_minus_spc}])

    return outdf


def recursive_threshold(v_channel, classes=5):
    try:
        threshs = threshold_multiotsu(v_channel, classes=classes)
        return threshs
    except ValueError as e:
        # if a ValueError occurs, print the error and try with one less class
        print(f"Error encountered: {e}. Trying with classes={classes-1}")
        if classes > 1:
            return recursive_threshold(v_channel, classes - 1)
        else:
            raise ValueError("threshold_multiotsu failed for all classes down to 1")


def threshold_img(img, classes=3):
    # get hsv value from image instead
    img_hsv = color.rgb2hsv(img)
    v_channel = img_hsv[:, :, 2]

    # if above fails for the default threshold of 5, recursively try 5 - n
    # for each increasing n + 1 until it works.
    threshs = recursive_threshold(v_channel, classes=classes)
    regions = np.digitize(v_channel, bins=threshs)
    # split out each category (use only the most intense-color one)
    maxval = np.unique(regions).max()
    region_max = deepcopy(regions) 
    region_max[np.where(region_max != maxval)] = 0

    # pick out first two categories
    categs = np.unique(regions)
    categs.sort()
    last2 = categs[-2:]
    region_use = deepcopy(regions) 
    idxs = np.where(np.isin(region_use, last2, invert=True))
    region_use[idxs] = 0

    return region_use


def main(to_process, savedir,
         close_radius, min_clean_size, min_watershed_dist, dapi_watershed, sigma, spc_watershed=None, spc_sigma=None, spc_dilate=None, spc_cleansize=None,
         dapi_hopx_watershed=None, dapi_aqp5_watershed=None, erdr1_sigma=None, dapi_colocalize=True, aqp5_sigma=None,
         dapi_img=None, hopx_img=None, tdt_img=None, spc_img=None, lcn2_img=None, krt8_img=None, aqp5_img=None, erdr1_img=None, gfp_img=None, rfp670_img=None, rab_img=None,
         plot=False, section=None, font=0.005, stopcheck=True):

    print("Will analyze " + to_process)
    os.makedirs(savedir, exist_ok=True)

    print("Parameters:")
    print("close radius: " + str(close_radius))
    print("min clean size: " + str(min_clean_size))
    print("general min watershed dist: " + str(min_watershed_dist))
    print("dapi hopx watershed dist: " + str(dapi_hopx_watershed))
    print("spc watershed: " + str(spc_watershed))
    print("gaussian blur sigma: " + str(sigma))
    print("spc dilate: " + str(spc_dilate))
    print("spc clean size: " + str(spc_cleansize))
    print("colocalizing with dapi: " + str(dapi_colocalize))


    all_imgs = {'dapi': dapi_img, 'hopx': hopx_img, 'tdt': tdt_img, 'spc': spc_img, 'lcn2': lcn2_img, 'krt8': krt8_img,
                'aqp5': aqp5_img, 'erdr1': erdr1_img, 'gfp': gfp_img, 'rfp670': rfp670_img, 'rab': rab_img}

    # check if any information in the image. if all blank, exit
    if stopcheck:
        if len(np.unique(all_imgs['dapi'])) == 1:
            raise ValueError("Image is empty--returning None")

    print("Thresholding images")
    thresh_imgs = {}
    for key, arr in all_imgs.items():
        if arr is not None:
            print("Thresholding " + key)
            thresh_imgs[key] = threshold_img(arr, classes=5)
        else:
            print(f"No image provided for {key}, skipping thresholding.")
            thresh_imgs[key] = None

    
    if to_process == "lcn2":
        print("Processing Lcn2")
        out = lcn2(lcn2_msk=thresh_imgs['lcn2'], tdt_msk=thresh_imgs['tdt'], spc_msk=thresh_imgs['spc'], dapi_msk=thresh_imgs['dapi'], plot=plot, savedir=savedir,
                         min_clean_size=min_clean_size, close_radius=close_radius, spc_dilate=spc_dilate, spc_cleansize=spc_cleansize,
                         spc_sigma=spc_sigma, sigma=sigma, spc_watershed=spc_watershed, dapi_watershed=dapi_watershed,
                         single_watershed=min_watershed_dist, section=section, font=font)

    elif to_process == "erdr1":
        print("Processing erdr1")
        out = erdr1(erdr1_msk=thresh_imgs['erdr1'], tdt_msk=thresh_imgs['tdt'], rfp670_msk=thresh_imgs['rfp670'], dapi_msk=thresh_imgs['dapi'], plot=plot, savedir=savedir,
                         min_clean_size=min_clean_size, close_radius=close_radius, spc_dilate=spc_dilate, erdr1_sigma=erdr1_sigma,
                         sigma=sigma, dapi_watershed=dapi_watershed, spc_watershed=spc_watershed, spc_cleansize=spc_cleansize, spc_sigma=spc_sigma, 
                         single_watershed=min_watershed_dist, section=section, font=font, dapi_colocalize=dapi_colocalize)
        
    elif to_process == "spc":
        print("Processing SPC")
        out = tdt(tdt_msk=thresh_imgs['tdt'], spc_msk=thresh_imgs['spc'], plot=plot, savedir=savedir,
                         close_radius=close_radius, min_clean_size=min_clean_size,
                            single_watershed=min_watershed_dist, spc_watershed=spc_watershed,
                            sigma=sigma, spc_sigma=spc_sigma, spc_dilate=spc_dilate, spc_cleansize=spc_cleansize,
                            section=section, font=font, run_externally=True, run_spc=True)
        
    elif to_process == "hopx":
        print("Processing Hopx")
        out = hopx(hopx_msk=thresh_imgs['hopx'], dapi_msk=thresh_imgs['dapi'], tdt_msk=thresh_imgs['tdt'], spc_msk=thresh_imgs['spc'], rab_msk=thresh_imgs['rab'], plot=plot, savedir=savedir,
                         sigma=sigma, close_radius=close_radius, min_clean_size=min_clean_size,
                         spc_cleansize=spc_cleansize, single_watershed=min_watershed_dist, dapi_watershed=dapi_watershed,
                         dapi_hopx_watershed=dapi_hopx_watershed, spc_watershed=spc_watershed, spc_dilate=spc_dilate, spc_sigma=spc_sigma,
                         dapi_sigma=sigma, hopx_spc_clean=min_clean_size, hopx_tdt_clean=min_clean_size,
                         section=section, font=font, dapi_colocalize=dapi_colocalize)
    elif to_process == "aqp5":
        print("Processing Aqp5")
        out = aqp5(aqp5_msk=thresh_imgs['aqp5'], dapi_msk=thresh_imgs['dapi'], tdt_msk=thresh_imgs['tdt'], spc_msk=thresh_imgs['spc'], plot=plot, savedir=savedir,
                         sigma=sigma, aqp5_sigma=aqp5_sigma, single_watershed=min_watershed_dist, dapi_watershed=dapi_watershed, dapi_aqp5_watershed=dapi_aqp5_watershed,
                         spc_watershed=spc_watershed, spc_sigma=spc_sigma, spc_cleansize=spc_cleansize, spc_dilate=spc_dilate,
                         close_radius=close_radius, min_clean_size=min_clean_size, aqp5_spc_clean=min_clean_size, aqp5_tdt_clean=min_clean_size, dapi_colocalize=dapi_colocalize,
                         section=section, font=font)

    else:
        raise ValueError("'to_process' must be 'tdt' 'lcn2', or 'erdr1'")

    return out
