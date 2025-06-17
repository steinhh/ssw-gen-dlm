; Produced by compile_sfit.pro                               
;                                                            
; To direct these automatically produced files to another    
; directory, set the environment variable IDL_COMPILE_DIR    
; to point at the directory.                                 
                                                             
PRO cf_p0_,x,a,f,pder                                        
    on_error,0                                               
                                                             
                                                             
                                                             
    nx = n_elements(x)                                       
                                                             
    use_pder = (n_params() eq 4)                             
                                                             
    type = datatype(a,2)                                     
    if use_pder then begin                                   
        pder = make_array(nx,n_elements(a),type=type,/nozero)
    end                                                      
                                                             
    f = make_array(nx,type=type)                             
    atemp = a(0:0)                                           
    if use_pder then begin                                   
       pder_temp = 1                                         
       comp_poly,x,atemp,ftemp,pder_temp                     
       pder(0,0) = pder_temp                                 
    end else begin                                           
       comp_poly,x,atemp,ftemp                               
    end                                                      
                                                             
    f = temporary(f) + ftemp                                 
end                                                          
