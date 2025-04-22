function DspRst (ite,Nsam,Nsze,Ncal)
if ite==1
    fprintf( '---------------------------------------------\n') ;
    fprintf( '     i      Nsam      Nc      Ns      Ncal   \n') ;
    fprintf( '---------------------------------------------\n') ;
end
fprintf( '%6d    %6d  %6d  %6d    %6d \n',ite,Nsam,Nsze(1),Nsze(2),Ncal);
return