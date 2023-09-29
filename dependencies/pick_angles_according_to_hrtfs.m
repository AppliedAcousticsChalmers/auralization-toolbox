function [azi_out_deg, ele_out_deg, indices_hrirs] = pick_angles_according_to_hrtfs(hrirs_sofa, azi_in_deg, ele_in_deg)
% returns the nearest neighbors of an angular grid that are comprised in
% the given HRTF set
%
% input: 
%    azi_in_deg, ele_in_deg: azimuth and elevtion in deg
%
% output: 
%    all in deg

indices_hrirs = SOFAfind(hrirs_sofa, azi_in_deg, ele_in_deg);

% quantize sdm_data to align with the HRTF angles
azi_out_deg = hrirs_sofa.SourcePosition(indices_hrirs, 1).';
ele_out_deg = hrirs_sofa.SourcePosition(indices_hrirs, 2).';

end