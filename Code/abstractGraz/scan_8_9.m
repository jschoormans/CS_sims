
clear all; close all; clc;
cd('/home/jschoormans/lood_storage/divi/Projects/cosart/scans/Graz1-phantoms_NSA_resolution')
files=dir('*.raw');

MR=MRecon(files(7).name)
MR.Parameter.Recon.ImmediateAveraging='yes'
MR.ReadData;
MR.RandomPhaseCorrection;
MR.PDACorrection;
MR.DcOffsetCorrection;
MR.MeasPhaseCorrection;
MR.SortData;
K=MR.Data{1};
MR.Perform
MR.ShowData



%%










