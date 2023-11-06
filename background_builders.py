#!/usr/bin/env python3

import numpy as np


class SingleBackgroundBuilder:
    #def __init__(self, bin_angle_masks, cv_img):
        #self.cv_img = cv_img
        #self.bin_angle_masks = bin_angle_masks

        # we'll use the background image to recompute the binning
        #self.bin_background = self.compute_bin_values(self.cv_img)
        #print("here1", self.bin_background)

    def __init__(self, bin_angle_masks, cv_img, bin_background):

        self.cv_img = cv_img
        self.bin_angle_masks = bin_angle_masks
        self.number_of_bins = len(self.bin_angle_masks["left"])
        #Left_Background = np.zeros((self.number_of_bins))
        Left_Background = bin_background[0:self.number_of_bins]
        #Right_Background = np.zeros((self.number_of_bins))
        Right_Background = bin_background[self.number_of_bins:self.number_of_bins*2]
        self.bin_background ={ "left": Left_Background, "right": Right_Background}

        #print("Left_Background", Left_Background)
        #print("Right_Background", Right_Background)

        #print("here1", self.number_of_bins)

    def compute_bin_values(self, img):
        raw_array_dict = {}
        for wing_key in ["left", "right"]:
            raw_array_list = []
            for _, single_bin_mask in self.bin_angle_masks[wing_key].items():
                raw_bin_average = img[single_bin_mask].mean()
                raw_array_list.append(raw_bin_average)
            raw_array_dict[wing_key] = np.array(raw_array_list)
        return raw_array_dict

    def reset_background(self, cv_img):
        self.cv_img = cv_img
        self.bin_background = self.compute_bin_values(cv_img)
        #print("here22", self.bin_background)

    def change_bins(self, bin_angle_masks, bin_background):
        self.bin_angle_masks = bin_angle_masks
        self.number_of_bins = len(self.bin_angle_masks["left"])
        #Left_Background = np.zeros((self.number_of_bins))
        Left_Background = bin_background[0:self.number_of_bins]
        #Right_Background = np.zeros((self.number_of_bins))
        Right_Background = bin_background[self.number_of_bins:self.number_of_bins*2]
        self.bin_background ={ "left": Left_Background, "right": Right_Background}

        #self.bin_angle_masks = bin_angle_masks
        #self.bin_background = self.compute_bin_values(self.cv_img)
        #print("here23", self.bin_background)
        #print("here12", len(self.bin_background["right"]))
        #print("here13", len(Right_Background))
        #print("here14", self.number_of_bins*2)

    def get_background_subtracted_array(self, new_img):
        bin_background_subtracted = {}
        raw_array_dict = self.compute_bin_values(new_img)
        for wing_key in ["left", "right"]:
            bin_background_subtracted[wing_key] = (
                raw_array_dict[wing_key] - self.bin_background[wing_key]
            )
        return bin_background_subtracted, self.bin_background


class MultipleBackgroundBuilder:
    def __init__(self, bin_angle_masks, cv_img, bin_background):
        self.bin_angle_masks = bin_angle_masks
        self.number_of_bins = len(self.bin_angle_masks["left"])
        self.row_idx = 0
        self.max_rows = 1000

        first_snapshot_background = self.compute_bin_values(cv_img)

        #first_snapshot_background_min = [min(first_snapshot_background)]    ################################# previous code
        #print(first_snapshot_background)

        self.raw_bin_backgrounds = {}
        for hinge_name in ["left", "right"]:
            self.raw_bin_backgrounds[hinge_name] = np.zeros((10, self.number_of_bins))
            self.raw_bin_backgrounds[hinge_name][:] = np.nan
            self.raw_bin_backgrounds[hinge_name][0, :] = first_snapshot_background[hinge_name]        # previous code

            #minimum = min(first_snapshot_background[hinge_name])        #################################
            #first_snapshot_background[hinge_name][:] = minimum        ################################# Use the minmum for each wing
            #print(minimum)              #####################################################
            
        self.counter = 0
        self.bin_background = first_snapshot_background
        #print('background',self.bin_background)              #####################################################

    def compute_bin_values(self, img):
        raw_array_dict = {}
        for wing_key in ["left", "right"]:
            raw_array_list = []
            for _, single_bin_mask in self.bin_angle_masks[wing_key].items():
                raw_bin_average = img[single_bin_mask].mean()
                raw_array_list.append(raw_bin_average)
            raw_array_dict[wing_key] = np.array(raw_array_list)
        return raw_array_dict

    def reset_background(self, cv_img):
        first_snapshot_background = self.compute_bin_values(cv_img)
        self.number_of_bins = len(self.bin_angle_masks["left"])

        self.raw_bin_backgrounds = {}
        for hinge_name in ["left", "right"]:
            self.raw_bin_backgrounds[hinge_name] = np.zeros((10, self.number_of_bins))
            self.raw_bin_backgrounds[hinge_name][:] = np.nan
            self.raw_bin_backgrounds[hinge_name][0, :] = first_snapshot_background[hinge_name]
        self.bin_background = first_snapshot_background

    def change_bins(self, bin_angle_masks, bin_background):
        self.bin_angle_masks = bin_angle_masks
        self.number_of_bins = len(self.bin_angle_masks["left"])

        # reset all backgrounds
        self.raw_bin_backgrounds = {}
        self.bin_background
        for hinge_name in ["left", "right"]:
            self.raw_bin_backgrounds[hinge_name] = np.zeros((10, self.number_of_bins))
            self.raw_bin_backgrounds[hinge_name][:] = np.nan
            self.bin_background[hinge_name] = np.zeros((self.number_of_bins))

    def get_background_subtracted_array(self, new_img):
        bin_backround_subtracted = {}

        raw_array_dict = self.compute_bin_values(new_img)
        #print("raw",raw_array_dict)   #####################################################################################

        if ((self.counter % 50) == 0) & (self.row_idx < self.max_rows):
            for wing_key in ["left", "right"]:
                self.row_idx = np.count_nonzero(
                    ~np.isnan(self.raw_bin_backgrounds[wing_key][:, 0])
                )
                if self.row_idx == len(self.raw_bin_backgrounds[wing_key]):
                    expand_mat = self.raw_bin_backgrounds[wing_key].copy()
                    expand_mat[:] = np.nan
                    self.raw_bin_backgrounds[wing_key] = np.concatenate(
                        [self.raw_bin_backgrounds[wing_key], expand_mat]
                    )

                self.raw_bin_backgrounds[wing_key][self.row_idx, :] = raw_array_dict[
                    wing_key
                ]
                self.bin_background[wing_key] = np.nanmin(
                    self.raw_bin_backgrounds[wing_key], axis=0
                )

        for wing_key in ["left", "right"]:
            bin_backround_subtracted[wing_key] = (
                raw_array_dict[wing_key] - self.bin_background[wing_key]
            )
        #print("raw",raw_array_dict)   #####################################################################################
        #print("back",self.bin_background)
        #print("diff",bin_backround_subtracted)

        self.counter += 1
        #print(self.bin_background)              #####################################################
        return bin_backround_subtracted, self.bin_background 
