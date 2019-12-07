function [symbol_out,bits]=QAM_slicer(OMEGA,BITS,slicer_input)

distance_table=abs(OMEGA-slicer_input*ones(size(OMEGA)));

[trash,min_index]=min(distance_table(:));

bits=BITS(min_index,:);
symbol_out=OMEGA(min_index,:);
