%EXP 4 JULY 2016: POMEGRANATE NOISY 2D SCANS WITH HIGH NSA
cd('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/Code')
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/MRIPhantomv0-8'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/tightfig'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/WaveLab850'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/sparseMRI_v0.2'))
%%
clear all; close all; clc;


folder='/home/jschoormans/lood_storage/divi/Ima/parrec/Jasper/pomegranate/'
file7='ph_04072016_1745477_7_2_wipt1wffeclearV4.raw' 
file12='ph_04072016_1759344_12_2_wipt1wffeclearV4.raw' 
MR=MRecon([folder,file7])
%%
MR.Parameter.Recon.ImmediateAveraging='No' 
MR.ReadData;
MR.DcOffsetCorrection;
MR.PDACorrection;
MR.RandomPhaseCorrection;
MR.MeasPhaseCorrection;
MR.SortData;
%%

