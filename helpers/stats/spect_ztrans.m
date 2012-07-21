function zval=spect_ztrans(data,beta,dof)
%z-transform for coherence measures
%
%

q=sqrt(-(dof-2).*log(1-data.^2));
zval=beta*(q-beta);

end


