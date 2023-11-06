#!/usr/bin/evn python3

import cv2
import numpy as np
from dataclasses import dataclass
import tkinter as tk
import time
import PyCapture2
import socket
#import selectors
import select

HOST = "127.0.0.1"  # Standard loopback interface address (localhost)
PORT = 65438  # Port to listen on (non-privileged ports are > 1023)

#Date = 't'
Date = input("Enter Date or 't' for Test or 'b' for new background: ")
if(Date == 't' or Date == 'T'):
    New_Background = 0
    Filename = 'Test'
    Real = 0
elif(Date == 'b' or Date == 'B'):
    New_Background = 1
    Filename = 'Test'
    Real = 0
else:
    Fly = input("Enter Fly Number: ")
    Trial = input("Enter Trial Number: ")
    Filename = Date + "_" + Fly + "_" + Trial
    Real = 1
print(Filename)

Input = 1               # 0 = video file, 1 = camera                     Take input from video file or camera
Debug = 0               # 0 = no debug, 1 = debug                        Save debug file
debug_BG = 0            # 0 = no debug, 1 = debug                        print background values
#Save = Real             # 0 = no save, 1 = save                          Save output file of time, left, right, frames
VideoOut = 1            # 0 = no video out, 1 = video out                Save video file
#New_Background = 0      # 0 = no new background, 1 = new background      Use 'Multiple' backgrounds to create a new background or use 'Single' to use saved background
Test = Real             # 0 = 'Test', 1 = enter new filename             Use "test" for the filename or enter new filename

if(Input):
    bus = PyCapture2.BusManager()
    numCams = bus.getNumOfCameras()
    camera = PyCapture2.Camera()
    uid = bus.getCameraFromIndex(0)
    camera.connect(uid)
    camera.startCapture()
    image = camera.retrieveBuffer() #Read out the first frame to check size..
    #row_bytes = float(len(image.getData())) / float(image.getRows());
    cv_image = np.array(image.getData(), dtype="uint8").reshape((image.getRows(), image.getCols()) );
    frame_width = cv_image.shape[1] #640  
    frame_height = cv_image.shape[0] #480 
    dim = (int(frame_width),int(frame_height)) #int is critical here
    fourcc = cv2.VideoWriter_fourcc(*'XVID')
    font                   = cv2.FONT_HERSHEY_COMPLEX
    bottomLeftCornerText = (25,450)
    fontScale              = 0.75 #was 1
    fontColor              = (255,255,255)
    lineType               = 2
    #Video_Folder = "Videos/"
    Video_Folder = "C:/Experiments/"
else:
    #Input_Video = "fc2_save_2023-05-19-1422"
    #Input_Video = "Video_0_0_0"
    Input_Video = "Video_1_4_g6"
    #Input_Video = "itzel_test_video"
    Video_Folder = "Videos/"
    #Video_Folder = "Videos/20230530/"
    #Location = Video_Folder + Input_Video + ".mp4"
    Location = Video_Folder + Input_Video + ".avi"
    video_capture = cv2.VideoCapture(Location, 0)
    total_n_frames = video_capture.get(cv2.CAP_PROP_FRAME_COUNT)
    print(total_n_frames)

def LoadFile(Location):         #************** load input array
    Input_Array = np.loadtxt(Location, delimiter =",")  
    return Input_Array

def FileSave(Output_File, Folder, Parameter_Array):
    np.savetxt(Folder + Output_File + '.csv', Parameter_Array, delimiter = ',')

#Param_Input = "Input_Parameters"
Param_Input = "Output_Parameters"
Param_Output = "Output_Parameters"
Back_File = "Background"
Debug_File = "Debug_" + Filename
Output_File = "Python_Output_" + Filename 
Param_File = "Python_Parameters_" + Filename
Background_File = "Background"
Param_Folder = "Parameters/"
Output_Folder = "C:/Experiments/"
Location_Input = Param_Folder + Param_Input + ".csv"
Location_Back = Param_Folder + Back_File + ".csv"
Input_Array = LoadFile(Location_Input)
Parameter_Array = np.empty([8])
Debug_Array3 = np.zeros([1,63]) ##number of bins +3
Save_Array = np.empty([5])
print(Input_Array)
Parameter_Array = Input_Array
Background_Array = LoadFile(Location_Back)
#print(Background_Array)

root = tk.Tk()
root.geometry("400x400")

label1 = tk.Label(root, text="X center:") 
label1.grid(row=0, column=0)
label2 = tk.Label(root, text="Y Center:") 
label2.grid(row=1, column=0)
label3 = tk.Label(root, text="Rotation:") 
label3.grid(row=2, column=0)
label4 = tk.Label(root, text="Torso Width:") 
label4.grid(row=3, column=0)
label5 = tk.Label(root, text="Inner Radius:") 
label5.grid(row=4, column=0)
label6 = tk.Label(root, text="Outer Radius:") 
label6.grid(row=5, column=0)
label7 = tk.Label(root, text="Bin Diff Thres:") 
label7.grid(row=6, column=0)
label8 = tk.Label(root, text="Flight Thres:") 
label8.grid(row=7, column=0)

slider_value1  = tk.DoubleVar()
slider_value2  = tk.DoubleVar()
slider_value3  = tk.DoubleVar()
slider_value4  = tk.DoubleVar()
slider_value5  = tk.DoubleVar()
slider_value6  = tk.DoubleVar()
slider_value7  = tk.DoubleVar()
slider_value8  = tk.DoubleVar()

slider1 = tk.Scale(root, variable=slider_value1 , from_=200, to=400, orient=tk.HORIZONTAL, length=300)  #variable=value1,
slider1.grid(row=0, column=1)
slider2 = tk.Scale(root, variable=slider_value2, from_=150, to=300, orient=tk.HORIZONTAL, length=300)
slider2.grid(row=1, column=1)
slider3 = tk.Scale(root, variable=slider_value3, from_=0, to=360, orient=tk.HORIZONTAL, length=300)
slider3.grid(row=2, column=1)
slider4 = tk.Scale(root, variable=slider_value4, from_=10, to=30, orient=tk.HORIZONTAL, length=300)
slider4.grid(row=3, column=1)
slider5 = tk.Scale(root, variable=slider_value5, from_=50, to=150, orient=tk.HORIZONTAL, length=300)
slider5.grid(row=4, column=1)
slider6 = tk.Scale(root, variable=slider_value6, from_=100, to=250, orient=tk.HORIZONTAL, length=300)
slider6.grid(row=5, column=1)
slider7 = tk.Scale(root, variable=slider_value7, from_=0, to=200, orient=tk.HORIZONTAL, length=300)
slider7.grid(row=6, column=1)
slider8 = tk.Scale(root, variable=slider_value8, from_=0, to=500, orient=tk.HORIZONTAL, length=300)
slider8.grid(row=7, column=1)

value1 = tk.IntVar()
value2 = tk.IntVar()
value3 = tk.IntVar()
value4 = tk.IntVar()
value5 = tk.IntVar()
value6 = tk.IntVar()
value7 = tk.IntVar()
value8 = tk.IntVar()

slider_value1.set(Input_Array[0]) 
slider_value2.set(Input_Array[1]) 
slider_value3.set(Input_Array[2]) 
slider_value4.set(Input_Array[3]) 
slider_value5.set(Input_Array[4]) 
slider_value6.set(Input_Array[5]) 
slider_value7.set(Input_Array[6]) 
slider_value8.set(Input_Array[7]) 

value1.set(Input_Array[0]) 
value2.set(Input_Array[1]) 
value3.set(Input_Array[2]) 
value4.set(Input_Array[3]) 
value5.set(Input_Array[4]) 
value6.set(Input_Array[5]) 
value7.set(Input_Array[6]) 
value8.set(Input_Array[7]) 

def update_value1(event):
    value1.set(slider1.get())
def update_value2(event):
    value2.set(slider2.get())
def update_value3(event):
    value3.set(slider3.get())
def update_value4(event):
    value4.set(slider4.get())
def update_value5(event):
    value5.set(slider5.get())
def update_value6(event):
    value6.set(slider6.get())
def update_value7(event):
    value7.set(slider7.get())
def update_value8(event):
    value8.set(slider8.get())

slider1.bind("<ButtonRelease-1>", update_value1)
slider2.bind("<ButtonRelease-1>", update_value2)
slider3.bind("<ButtonRelease-1>", update_value3)
slider4.bind("<ButtonRelease-1>", update_value4)
slider5.bind("<ButtonRelease-1>", update_value5)
slider6.bind("<ButtonRelease-1>", update_value6)
slider7.bind("<ButtonRelease-1>", update_value7)
slider8.bind("<ButtonRelease-1>", update_value8)

from background_builders import (SingleBackgroundBuilder, MultipleBackgroundBuilder,)

DEG2RAD = np.pi / 180
RAD2DEG = 180 / np.pi

@dataclass
class EdgeParameters:
    bin_diff_threshold: float = 0.1
    flight_no_flight_threshold: float = 1.0
    bottom_view: bool = False

@dataclass
class ViewParameters:
    show_mask: bool = True
    show_bin_values: bool = True

@dataclass
class BackgroundParameters:
    background_builder_type: str = "multiple"
    #background_builder_type: str = "single"
    reset_background: bool = False

@dataclass
class BinningParameters:
    inner_radius: int
    outer_radius: int
    center_x: int
    center_y: int
    torso_width: int = 20
    number_of_bins: int = 20
    min_angle: int = 60
    max_angle: int = 180
    body_angle: int = 0

def find_threshold_crossing(array, threshold):
    threshold_crossing = 0
    dx = 1
    for idx, value in enumerate(array[1:], start=1):
        prior_value = array[idx - 1]
        if (prior_value < threshold):
            threshold_crossing = 0
            break
        elif ((prior_value > threshold) & (value < threshold)):
            dy = value - prior_value
            dydx = dy / dx
            threshold_crossing = (idx - 1) + (threshold - prior_value) / dydx
            break
            
    #if(threshold_crossing >1):
    #    print('crossings', threshold, array, threshold_crossing)      #########################################
    return threshold_crossing

def find_threshold_crossing2(array, threshold):
    threshold_crossings = []
    dx = 1
    for idx, value in enumerate(array[1:], start=1):
        prior_value = array[idx - 1]

        if ((prior_value < threshold) & (value > threshold)) | (
            (prior_value > threshold) & (value < threshold)
        ):
            dy = value - prior_value
            dydx = dy / dx
            xt = (idx - 1) + (threshold - prior_value) / dydx
            threshold_crossings.append(xt)
    #print('crossings', threshold, array, threshold_crossings)      #########################################
    return np.array(threshold_crossings)

class LeadingEdgeFinder:
    def __init__(
        self,
        bin_params=None,
        edge_params=EdgeParameters(),
        background_params=BackgroundParameters(),
        view_params=ViewParameters(),
    ):
        self._edge_params = edge_params
        self._bin_params = bin_params
        self._background_params = background_params

        self._view_params = view_params
        self._bin_edges = None

        self._display_masks = None
        self._reset_binning_flag = False
        self._reset_background_flag = False
        self._latest_debug_data = None

    def change_params(self, param_group):

        param_type = type(param_group)

        if param_type == BinningParameters:
            # reset bins computation in upcoming image
            self._bin_params = param_group
            self._reset_binning_flag = True

        elif param_type == BackgroundParameters:
            # reset background next image
            self._background_params = param_group
            self._reset_background_flag = True

        elif param_type == ViewParameters:
            self._view_params = param_group

        elif param_type == EdgeParameters:
            self._edge_params = param_group
        else:
            raise Exception("Not a recognized parameter to set")

    def _first_frame_init(self, cv_img):

        self._display_masks = [np.zeros_like(cv_img), np.zeros_like(cv_img)]
        self.img_height, self.img_width = cv_img.shape[:2]

        if not self._bin_params:
            min_image_dim = min(self.img_height, self.img_width)
            half_len = min_image_dim // 2
            third_len = min_image_dim // 3
            quarter_len = min_image_dim // 4

            self._bin_params = BinningParameters(
                quarter_len, third_len, half_len, half_len
            )
            self._bin_angle_masks = self._compute_bin_masks()

        self._display_masks = [np.zeros_like(cv_img), np.zeros_like(cv_img)]
        self._bin_angle_masks = self._compute_bin_masks()
        self._initialize_background_builder(cv_img,Background_Array)
        self.background_builder.change_bins(self._bin_angle_masks, Background_Array)

    def _initialize_background_builder(self, cv_img,Background_Array):

        if self._background_params.background_builder_type == "single":
            self.background_builder = SingleBackgroundBuilder(
                self._bin_angle_masks, cv_img, Background_Array
            )
        elif self._background_params.background_builder_type == "multiple":
            self.background_builder = MultipleBackgroundBuilder(
                self._bin_angle_masks, cv_img, Background_Array
            )

    def analyze_image(self, cv_img):
        Debug_Array = np.empty([42])

        if self._display_masks is None:
            self._first_frame_init(cv_img)

        if self._reset_binning_flag:
            self._bin_angle_masks = self._compute_bin_masks()
            self.background_builder.change_bins(self._bin_angle_masks, Background_Array)
            self._reset_binning_flag = False

        if self._reset_background_flag:
            self._initialize_background_builder(cv_img)
            self.background_builder.reset_background(cv_img)
            self._reset_background_flag = False

        # do background subtraction
        # We only subtract the mean background per bin from the mean bin intensity
        # so we don't have to diff any pixels not in our bin mask
        bins_background_sub, background = self.background_builder.get_background_subtracted_array(cv_img)
        #Left_Bins = np.array(bins_background_sub['left'])   
        #Right_Bins = np.array(bins_background_sub['right'])
        Left_Bins = np.array(background['left'])   
        Right_Bins = np.array(background['right'])
        array_max = (np.mean(Left_Bins) + np.mean(Right_Bins))/2
        if(debug_BG):
            Left_bg_min = str(round(np.min(Left_Bins),2))
            Right_bg_min = str(round(np.min(Right_Bins),2))
            Left_bg_max = str(round(np.max(Left_Bins),2))
            Right_bg_max = str(round(np.max(Right_Bins),2))
            Left_sub_min = str(round(np.min(bins_background_sub['left']),2))
            Right_sub_min = str(round(np.min(bins_background_sub['right']),2))
            print("bg_min: ", Left_bg_min, Right_bg_min, "    bg_max", Left_bg_max, Right_bg_max, "    sub_min: ", Left_sub_min, Right_sub_min)
        #print("Mean",array_max)  #*************************************************************************************************************

        # Find average bin intensity and normalized bin intensity (for flight/no flight and leading edge position respectively)
        avg_bin_intensity = {}
        norm_bin_intensity = {}
        for wing_key in ["left", "right"]:
            avg_bin_intensity[wing_key] = (
                bins_background_sub[wing_key].sum() / self._bin_params.number_of_bins
            )
            positive_bin_intensity = bins_background_sub[wing_key] - np.nanmin(
                bins_background_sub[wing_key]
            )
            #array_max = np.nanmax(positive_bin_intensity)      
            #print(wing_key, array_max, bins_background_sub, avg_bin_intensity[wing_key],self._edge_params.flight_no_flight_threshold)  #######################
            #print('minimum',np.nanmin(bins_background_sub[wing_key]))  #######################
            #array_max = background[wing_key])
            if array_max == 0:
                array_max = 1
            
            positive_bin_intensity = positive_bin_intensity / array_max
            norm_bin_intensity[wing_key] = positive_bin_intensity
            #print("help",wing_key,array_max, positive_bin_intensity)  ############################################################################

        # find bin index where threshold has been crossed
        # for wing angle, we will interpolate to get sub-bin accuracy
        idx_angle_bin = {}
        wing_angles = {}
        for wing_key in ["left", "right"]:
            #diff_idx = -np.diff(norm_bin_intensity[wing_key])

            threshold_crossing = find_threshold_crossing(
                norm_bin_intensity[wing_key], self._edge_params.bin_diff_threshold
            )
            """
            if len(threshold_array) > 0:
                threshold_sub_idx = threshold_array[-1]
            else:
                threshold_sub_idx = 0
            """

            # interpolate between bins to find sub bin precision pixels
            base_idx = int(threshold_crossing)
            leftover_idx = threshold_crossing - base_idx
            angle_float = (
                self._bin_edges[base_idx] + self.bin_angle_delta * leftover_idx
            )

            # store index of angle to plot in debug window
            idx_angle_bin[wing_key] = base_idx

            # store angle
            wing_angles[wing_key] = angle_float
        self._latest_debug_data = {
            "img": cv_img,
            "wing_angles": wing_angles,
            "avg_bin_intensity": avg_bin_intensity,
            "norm_bin_intensity": norm_bin_intensity,
            "idx_angle_bin": idx_angle_bin,
            #"Left_bg_min": Left_bg_min,
            #"Right_bg_min": Right_bg_min,
            #"Left_sub_min": Left_sub_min,
            #"Right_sub_min": Right_sub_min,
        }

        # If both values are below the threshold, set wing angles to zero, meaning no flight.
        if (
            avg_bin_intensity["left"] < self._edge_params.flight_no_flight_threshold
        ) and (
            avg_bin_intensity["right"] < self._edge_params.flight_no_flight_threshold
        ):
            wing_angles["left"] = 0.0
            wing_angles["right"] = 0.0

        # we assume view from bottom
        # if we have a top view, swap the angles
        if self._edge_params.bottom_view == False:
            wing_angles = {"right": wing_angles["left"], "left": wing_angles["right"]}
            avg_bin_intensity = {
                "right": avg_bin_intensity["left"],
                "left": avg_bin_intensity["right"],
            }
        #print('both',Left_Bins, Right_Bins)
        result = np.concatenate((Left_Bins, Right_Bins), axis=0)    #wing_angles["left"],Right_Bins,wing_angles["right"]
        result2 = np.array([wing_angles["left"] , wing_angles["right"]])
        #result2 = np.concatenate((result, wing_angles["left"] , wing_angles["right"]))
        #Debug_Array = np.concatenate((result, result2), axis=0)
        #print('conc', len(result),result, wing_angles["left"] , wing_angles["right"])
        #print('debug',len(Debug_Array),Debug_Array)

        return wing_angles, avg_bin_intensity, result, result2

    def _compute_bin_masks(self):

        width_idx_matrix, height_idx_matrix = np.meshgrid(
            range(self.img_width), range(self.img_height)
        )

        if self._bin_params.number_of_bins < 3:
            self._reset_binning_flag = False
            return

        self._bin_edges = np.linspace(
            self._bin_params.min_angle,
            self._bin_params.max_angle,
            self._bin_params.number_of_bins + 1,
        )

        self.bin_angle_delta = self._bin_edges[1] - self._bin_edges[0]
        self._bin_edges_dict = {}

        body_angle_rad = self._bin_params.body_angle / 180 * np.pi
        hingepoints_offset = 3 * np.pi / 4
        left_hingepoint = {}
        left_hingepoint["x"] = self._bin_params.center_x + int(
            self._bin_params.torso_width
            * (
                np.cos(body_angle_rad + hingepoints_offset)
                - np.sin(body_angle_rad + hingepoints_offset)
            )
        )
        left_hingepoint["y"] = self._bin_params.center_y + int(
            self._bin_params.torso_width
            * (
                np.sin(body_angle_rad + hingepoints_offset)
                + np.cos(body_angle_rad + hingepoints_offset)
            )
        )

        right_hingepoint = {}
        right_hingepoint["x"] = self._bin_params.center_x - int(
            self._bin_params.torso_width
            * (
                np.cos(body_angle_rad + hingepoints_offset)
                - np.sin(body_angle_rad + hingepoints_offset)
            )
        )

        right_hingepoint["y"] = self._bin_params.center_y - int(
            self._bin_params.torso_width
            * (
                np.sin(body_angle_rad + hingepoints_offset)
                + np.cos(body_angle_rad + hingepoints_offset)
            )
        )

        self.hinges = {"left": left_hingepoint, "right": right_hingepoint}

        # create empty mask
        bin_angle_masks = {}
        for wing_key in ["left", "right"]:

            # center matrices
            idx_width = width_idx_matrix - self.hinges[wing_key]["x"]
            idx_height = height_idx_matrix - self.hinges[wing_key]["y"]

            # compute radial and angular position of every pixel in the image, w.r.t. each wing
            # angle is zero at the tail, 180 at the head
            # the angle goes counter  clockwise for the left wing, clockwise for the right
            radial_mat = np.sqrt(idx_width ** 2 + idx_height ** 2)
            if wing_key == "left":
                angular_mat = (
                    np.arctan2(-idx_height, idx_width)
                    - np.pi / 2
                    + self._bin_params.body_angle * DEG2RAD
                )
                angular_mat = np.mod(angular_mat + np.pi, 2 * np.pi) - np.pi

            elif wing_key == "right":
                angular_mat = (
                    np.arctan2(-idx_height, -idx_width)
                    - np.pi / 2
                    - self._bin_params.body_angle * DEG2RAD
                )
                angular_mat = np.mod(angular_mat + np.pi, 2 * np.pi) - np.pi

            bin_angle_masks[wing_key] = {}
            for bin_idx in range(len(self._bin_edges) - 1):
                angle_rad_start = self._bin_edges[bin_idx] * DEG2RAD
                angle_rad_end = self._bin_edges[bin_idx + 1] * DEG2RAD
                in_slice = (
                    (radial_mat > self._bin_params.inner_radius)
                    & (radial_mat < self._bin_params.outer_radius)
                    & (angular_mat >= angle_rad_start)
                    & (angular_mat < angle_rad_end)
                )
                bin_angle_masks[wing_key][bin_idx] = in_slice

        # set display masks to zero prior to adding values
        self._display_masks[0] *= 0
        self._display_masks[1] *= 0

        counter = 0
        for wing_key in ["left", "right"]:
            for _, mask in bin_angle_masks[wing_key].items():
                if counter % 2 == 0:
                    self._display_masks[0][mask] = 20
                else:
                    self._display_masks[1][mask] = 10
                counter += 1

        return bin_angle_masks

    def create_debug_image(self):
        if not self._latest_debug_data:
            return

        angles = self._latest_debug_data["wing_angles"]
        normalized_array_dict = self._latest_debug_data["norm_bin_intensity"]
        cv_img = self._latest_debug_data["img"]
        idx_angle_bin = self._latest_debug_data["idx_angle_bin"]
        #Left_bg_min = self._latest_debug_data["Left_bg_min"]
        #Right_bg_min = self._latest_debug_data["Right_bg_min"]
        #Left_bg_max = self._latest_debug_data["Left_bg_max"]
        #Right_bg_max = self._latest_debug_data["Right_bg_max"]
        #Left_sub_min = self._latest_debug_data["Left_sub_min"]
        #Right_sub_min = self._latest_debug_data["Right_sub_min"]

        for wing_key in angles:
            if wing_key == "left":
                outmost_point = (
                    int(
                        self.hinges[wing_key]["x"]
                        + self._bin_params.outer_radius
                        * np.sin(
                            np.pi
                            + (angles[wing_key] - self._bin_params.body_angle) * DEG2RAD
                        )
                    ),
                    int(
                        self.hinges[wing_key]["y"]
                        + self._bin_params.outer_radius
                        * np.cos(
                            np.pi
                            + (angles[wing_key] - self._bin_params.body_angle) * DEG2RAD
                        )
                    ),
                )
            elif wing_key == "right":
                outmost_point = (
                    int(
                        self.hinges[wing_key]["x"]
                        + self._bin_params.outer_radius
                        * np.sin(
                            np.pi
                            - (angles[wing_key] + self._bin_params.body_angle) * DEG2RAD
                        )
                    ),
                    int(
                        self.hinges[wing_key]["y"]
                        + self._bin_params.outer_radius
                        * np.cos(
                            np.pi
                            - (angles[wing_key] + self._bin_params.body_angle) * DEG2RAD
                        )
                    ),
                )

            cv2.line(
                cv_img,
                (self.hinges[wing_key]["x"], self.hinges[wing_key]["y"]),
                outmost_point,
                (200, 200, 200),
                thickness=2,
            )

        for wing_key in angles:
            if self._view_params.show_bin_values:
                for bin_idx, bin_value in enumerate(normalized_array_dict[wing_key]):
                    if np.isnan(bin_value):
                        bin_value = 0
                    else:
                        bin_value = bin_value

                    y_position = int(50 + 10 * bin_idx)
                    if wing_key == "right":
                        x_0 = 20
                        x_value = int(x_0 + bin_value * self.img_width / 15)
                    elif wing_key == "left":
                        x_0 = self.img_width - 50
                        x_value = int(x_0 - bin_value * self.img_width / 15)
                    if bin_idx == idx_angle_bin[wing_key]:
                        color = (255, 255, 255)
                    else:
                        color = (150, 150, 150)
                    cv2.line(
                        cv_img,
                        (x_0, y_position),
                        (x_value, y_position),
                        color,
                        thickness=8,
                    )
        
        #print("bg_min: ", Left_bg_min, Right_bg_min, "bg_max: ", Left_bg_max, Right_bg_max, "sub_min: ", Left_sub_min, Right_sub_min)

        if self._view_params.show_mask:
            cv_img += self._display_masks[0]
            cv_img += self._display_masks[1]

        return cv_img

#def main(): ###############################################################################################
    """This is to test leading_edge_finder without any ROS packages.
    You will likely have to change the package path on your computer
    """

bin_params = BinningParameters(
    70,
    150,
    318,
    225,
    torso_width=25,
    number_of_bins=30,
    min_angle=60,
    max_angle=180,
    body_angle=250,
)

edge_params = EdgeParameters(
    bin_diff_threshold=0.1,
    flight_no_flight_threshold=1,
    bottom_view=False,
)

if(New_Background):
    background_params = BackgroundParameters(background_builder_type="multiple", reset_background=False)      ########################## set background type 
else:
    background_params = BackgroundParameters(background_builder_type="single", reset_background=False)
view_params = ViewParameters(show_bin_values=True, show_mask=True)
edge_finder = LeadingEdgeFinder(
    bin_params=bin_params,
    edge_params=edge_params,
    background_params=background_params,
    view_params=view_params,
)

Lastmarktime = time.time()
frameNo = 0
MaxLooptime = 0
DeltaTime = 0
LastPrintTime = time.time()
LastPrintFrame = 0
Cum_Wing_Diff = 0
Status = 0
LastUDPtime = time.time()
LastUDPCheck = time.time()

with socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as sTx:
    sRx = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sRx.bind((HOST, PORT+1))
    #sRx.setblocking(0)
    #sRx.settimeout(0.001)
    #sel = selectors.DefaultSelector()
    #sel.register(sRx, selectors.EVENT_READ, data=None)



    while True:
        Marktime = time.time()                  #check for time round the loop
        Looptime = Marktime - Lastmarktime
        Lastmarktime = Marktime
        if Looptime > MaxLooptime:
            MaxLooptime = Looptime

        if(Status == 1):
            if(round(Marktime - LastUDPCheck,1) > 3):
                readable, _, _ = select.select([sRx], [], [], 0.001) 
                LastUDPCheck = Marktime
                print('UDP Check')
        else:
            readable, _, _ = select.select([sRx], [], [], 0.001) 
            LastUDPCheck = Marktime
            #print('UDP Checking')
        #readable = sel.select(timeout=0.001)
        #print(time.time() - Marktime)      ######################################### time to check UDP
        if readable:
            data, addr = sRx.recvfrom(1024)
            print(data.decode())
            while(True):
                readable, _, _ = select.select([sRx], [], [], 0.001)
                if readable:
                    data, addr = sRx.recvfrom(1024)
                    print(data.decode())
                else:
                    break    
            if(str(data) == "b'Stop'"):##################### Stop Trial #################
                Status = 0
                print('Trial Stopped')
                while(True):
                    readable, _, _ = select.select([sRx], [], [], 0.001)
                    if readable:
                        data, addr = sRx.recvfrom(1024)
                        print(data.decode())
                    else:
                        break
                FileSave(Param_File, Output_Folder, Parameter_Array)
                if(1):
                    FileSave(Output_File, Output_Folder, Save_Array) 
                if(Input):
                    print(len(image_data_list))
                    for frame in image_data_list:
                        out.write(cv_image)
                    out.release()
            elif(Status):       ###################################   Trial Continues #################
                LastUDPtime = time.time()
                readable = 0
                print('Trial Continues')
            else:       ###################################   initialise Arrays - Start Trial #################
                print('Trial Started')
                Status = 1
                LastUDPtime = time.time()
                Filename = data.decode()
                frameNo = 0
                readable = 0
                Param_File = "Python_Par_" + Filename
                Output_File = "Python_Data_" + Filename
                Videofilename = 'Video_' + Filename
                if(Input):
                    image_data_list = []
                    out = cv2.VideoWriter(Video_Folder + Videofilename + ".avi",fourcc,30.0,dim,0)
                print(time.time() - Marktime)
        
        if(True):                               #check for GUI changes
            root.update()
            if (bin_params.center_x != value1.get()):
                bin_params.center_x = value1.get()
                Parameter_Array[0] = value1.get()
                edge_finder.change_params(bin_params)
            
            if (bin_params.center_y != value2.get()):
                bin_params.center_y = value2.get()
                Parameter_Array[1] = value2.get()
                edge_finder.change_params(bin_params)

            if (bin_params.body_angle != value3.get()):
                bin_params.body_angle = value3.get()
                Parameter_Array[2] = value3.get()
                edge_finder.change_params(bin_params)
            
            if (bin_params.torso_width != value4.get()):
                bin_params.torso_width = value4.get()
                Parameter_Array[3] = value4.get()
                edge_finder.change_params(bin_params)

            if (bin_params.inner_radius != value5.get()):
                bin_params.inner_radius = value5.get()
                Parameter_Array[4] = value5.get()
                edge_finder.change_params(bin_params)
            
            if (bin_params.outer_radius != value6.get()):
                bin_params.outer_radius = value6.get()
                Parameter_Array[5] = value6.get()
                edge_finder.change_params(bin_params)

            if (edge_params.bin_diff_threshold != value7.get()/1000):
                edge_params.bin_diff_threshold = value7.get()/1000
                Parameter_Array[6] = value7.get()
                edge_finder.change_params(edge_params)
            
            if (edge_params.flight_no_flight_threshold != value8.get()/100):
                edge_params.flight_no_flight_threshold = value8.get()/100
                Parameter_Array[7] = value8.get()
                edge_finder.change_params(edge_params)
        
        frameNo = 1 + frameNo
        if(Input):
            image = camera.retrieveBuffer()
            cv_image = np.array(image.getData(), dtype="uint8").reshape((image.getRows(), image.getCols()) );
            gray_frame = cv_image
            if(VideoOut):
                #cv2.putText(cv_image, 'Frames =' + str(frameNo), bottomLeftCornerText, font, fontScale, fontColor, lineType)
                if(Status):
                    txtlen = int(round(len(str(abs(frameNo)))-4)*20)
                    bottomLeftCornerOfText_3 = (int(frame_width-250-txtlen),int(50))
                    cv2.putText(cv_image, 'Frames =' + str(frameNo), bottomLeftCornerText, font, fontScale, fontColor, lineType)
                    image_data_list.append(cv_image)
                    #print('frame',frameNo)
        else:
            ret, frame = video_capture.read()
            gray_frame = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
            if(False):
                if video_capture.get(cv2.CAP_PROP_POS_FRAMES) == total_n_frames:
                    video_capture.set(cv2.CAP_PROP_POS_FRAMES, 0) ########################################################## restart video
            else:   
                if video_capture.get(cv2.CAP_PROP_POS_FRAMES) == 1800:
                    video_capture.set(cv2.CAP_PROP_POS_FRAMES, 1000) ########################################################## loop through these frames
                if video_capture.get(cv2.CAP_PROP_POS_FRAMES) == 1: 
                    video_capture.set(cv2.CAP_PROP_POS_FRAMES, 1000) ########################################################## start advance
                #print(video_capture.get(cv2.CAP_PROP_POS_FRAMES))
                #time.sleep(2)     ######################################## slow down

        angles, avg_bin_intensity, Background_Array, Angle_Array = edge_finder.analyze_image(gray_frame)
        Sendtime = round(time.time()*1000)
        aSendtime = np.array([Sendtime])
        sTx.sendto(str(Sendtime).encode(), (HOST, PORT))
        #Time = np.array([Marktime])
        Frame = np.array([frameNo])
        sFrame = str(Frame)
        Wing_Diff = np.array([Angle_Array[0] - Angle_Array[1]])
        #print(Angle_Array[0],Angle_Array[1])
        if(Angle_Array[0] ==0 or Angle_Array[1]==0):
            Wing_Diff = np.array([Angle_Array[0] - Angle_Array[1]+400]) #################### if either wing is not flying send 400 as dWLE - decode at the other end
        sWing_Diff = str(round(np.min(Wing_Diff),0))
        sTx.sendto(sWing_Diff.encode(), (HOST, PORT))
        Cum_Wing_Diff += Wing_Diff
        
        if(Frame%30 == 1):
            DeltaFrames = frameNo - LastPrintFrame
            if DeltaFrames < 1:
                DeltaFrames = 1
            LastPrintFrame = frameNo
            DeltaTime = Marktime - LastPrintTime
            Av_Wing_Diff = round(np.min(Cum_Wing_Diff/DeltaFrames),0)
            print('Wing Diff', Av_Wing_Diff,'  Frame', Frame,' Loop max:', round(MaxLooptime*1000), ' Av:' , round(DeltaTime*1000/DeltaFrames,1))
            #print('Status', Status, '  UDPcheck time', round(Marktime - LastUDPCheck,3))
            MaxLooptime = 0
            LastPrintTime = Marktime
            Cum_Wing_Diff = 0
    
        debug_image = edge_finder.create_debug_image()
        cv2.imshow(f"debug image", debug_image)
        
        if(Debug):
            Debug_Array = np.concatenate((Background_Array, Angle_Array), axis=0)
            Time = np.array([Marktime])
            Debug_Array2 = np.concatenate((Time, Debug_Array), axis=0)
            if(Debug_Array3.all == 0):
                Debug_Array3 = Debug_Array2
            else:
                Debug_Array3 = np.vstack((Debug_Array3, Debug_Array2))

        Debug_Array4 = np.concatenate((aSendtime, Wing_Diff, Angle_Array, Frame), axis=0)
        if(Status):
            if(frameNo == 1):
                Save_Array = Debug_Array4
            else:
                Save_Array = np.vstack((Save_Array, Debug_Array4))
         
        if(((LastUDPtime + 10) < Marktime)and(Status == 1)):
            Status = 0
            print('Trial Finished')
            while(True):
                    readable, _, _ = select.select([sRx], [], [], 0.001)
                    if readable:
                        data, addr = sRx.recvfrom(1024)
                        print(data.decode())
                    else:
                        break
            FileSave(Param_File, Output_Folder, Parameter_Array)
            if(1):
                FileSave(Output_File, Output_Folder, Save_Array) 
        if cv2.waitKey(1) & 0xFF == ord('q'):
        #if (True):# cv2.waitKey(1) & 0xFF == ord('q'):
            FileSave(Param_Output, Param_Folder, Parameter_Array)
            #FileSave(Param_File, Output_Folder, Parameter_Array)
            if(Debug):
                FileSave(Debug_File, Param_Folder, Debug_Array3)  
            #if(Save):
            #    FileSave(Output_File, Output_Folder, Save_Array)  
            if(New_Background):
                FileSave(Background_File, Param_Folder, Background_Array)
            if(Input):
                camera.stopCapture()
                camera.disconnect()
                #tout.release()
            break
sTx.close()
sRx.close()
cv2.destroyAllWindows()