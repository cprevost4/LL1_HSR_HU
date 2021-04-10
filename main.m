%        Hyperspectral super-resolution with variable images:             %
%              LL1-based recovery and blind unmixing                      %
%-------------------------------------------------------------------------%

% Copyright (c) 2021 Clemence Prevost, Ricardo A. Borsoi,
% Konstantin Usevich, David Brie, Jose ?C. M. Bermudez, Cedric Richard
% https://github.com/cprevost4/CCRB_Software
% Contact: clemence.prevost@univ-lorraine.fr

% This software reproduces the results from the paper called:
% "Hyperspectral super-resolution with variable images: LL1-based recovery 
% and blind unmixing" - C.Prévost, R. A. Borsoi, K.Usevich,
% D.Brie, J. M. Bermudez, C. Richard.
%
% In order to run the demo, you will need to add to your MATLAB path:
% - Tensorlab 3.0: https://www.tensorlab.net
%
%-------------------------------------------------------------------------%
%                              CONTENT
% - /baseline_algorithms: contains codes and adaptors for other methods
% - /data : contains data for synthetic examples (Section VI.D)
% - /demos : contains demo files that produce tables and figures
% - /figures : where the figures are saved
% - /metrics : contains the metrics use for comparison in the paper
% - /src : contains helpful files to run the demos
%
%-------------------------------------------------------------------------%
%                                MENU
% You can launch a specified demo by typing its number. The resulting tables
% and figures produced will be stored in the figures folder.
%
% 1:  produces Fig. 1 and Tables I, II 
% 2:  produces Fig. 2 and Tables III, IV 
% 3:  produces Fig. 3 and Table V 
% 4:  produces Fig. 4 and 5
% 5:  produces Fig. 6 and 7
% 6:  produces Fig. 8 and 9
% 7:  produces Fig. 11 and 12 
%
%-------------------------------------------------------------------------%

list_demos = ["fusion_semireal1" "fusion_semireal2" "fusion_semireal_novar"
     "unmixing_synthetic1" "unmixing_synthetic2" "unmixing_semireal1"
     "unmixing_semireal2"
    ];

prompt = "Which file do you want to run ?";
num = input(prompt);
eval(list_demos(num));