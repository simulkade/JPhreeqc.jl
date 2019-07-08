"""
function simple_seawater(;id::Int=1, temp_C::Real=25.0)
returns a simple phreeqc string that defines seawater composition 
In the original definition, charge keyword was not present. I added
it for the sake of charge balance
"""
function simple_seawater(;id::Int=1, temp_C::Real=25.0)
    return """
    SOLUTION $id    
        units   ppm
        pH      8.22  charge
        pe      8.451
        density 1.023
        temp    $temp_C
        Ca              412.3
        Mg              1291.8
        Na              10768.0
        K               399.1
        Si              4.28
        Cl              19353.0
        Alkalinity      141.682 as HCO3
        S(6)            2712.0
    END
    """
end # simple_seawater

function simple_brine(;id::Int=1, temp_C::Real=25.0)
    return """
    SOLUTION $id
        units mol/khw
        pH  7.0 charge
        Na  0.5
        Cl  0.5
    END
    """
end #simple_brine

function cd_music_calcite(;log_k=[-3.58, -2.8, -2.6, 12.85, -24.73, 10.15, 1.55, 0.96],
    dh=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    return """
    SURFACE_MASTER_SPECIES
        Chalk_a Chalk_aOH-0.667
        Chalk_c Chalk_cH+0.667
    SURFACE_SPECIES
        Chalk_cH+0.667 = Chalk_cH+0.667
        log_k 0.0
        -cd_music 0 0 0 0 0
        Chalk_aOH-0.667 = Chalk_aOH-0.667
        log_k 0.0
        -cd_music 0 0 0 0 0
        Chalk_cH+0.667 = Chalk_c-0.333 + H+
        log_k $(log_k[1]) # -3.58
        delta_h      $(dh[1])
        -cd_music -1 0 0 0 0
        Chalk_cH+0.667 + Ca+2 = Chalk_cCa+1.667 + H+
        log_k $(log_k[2]) # -2.8
        delta_h      $(dh[2])
        -cd_music -1 2 0 0 0
        Chalk_cH+0.667 + Mg+2 = Chalk_cMg+1.667 + H+
        log_k $(log_k[3]) # -2.6
        delta_h      $(dh[3])
        -cd_music -1 2 0 0 0
        Chalk_aOH-0.667 + H+ = Chalk_aOH2+0.333
        log_k $(log_k[4]) # 12.85
        delta_h      $(dh[4])
        -cd_music 1 0 0 0 0
        Chalk_aOH-0.667 = Chalk_aO-1.667 + H+
        log_k $(log_k[5]) # -24.73
        delta_h      $(dh[5])
        -cd_music -1 0 0 0 0
        Chalk_aOH-0.667 + CO3-2 + H+ = Chalk_aHCO3-0.667 + OH-
        log_k $(log_k[6]) # 10.15
        delta_h      $(dh[6])
        -cd_music .6 -.6 0 0 0
        Chalk_aOH-0.667 + CO3-2 = Chalk_aCO3-1.667 + OH-
        log_k $(log_k[7]) # 1.55
        delta_h      $(dh[7])
        -cd_music .6 -1.6 0 0 0
        Chalk_aOH-0.667 + SO4-2 = Chalk_aSO4-1.667 + OH-
        log_k $(log_k[8]) # 1.55
        delta_h      $(dh[8])
        -cd_music .6 -1.6 0 0 0
        # 2Chalk_aOH-0.667 + CaPO4- = Chalk_a2PO4Ca-0.334 + 2OH-
        # log_k 4.00
        # -cd_music  0 -4 0 .24 5
        # 2Chalk_aOH-0.667 + CaHPO4 = Chalk_a2HPO4Ca+0.666 + 2OH-
        # log_k 2.69
        # -cd_music  0 -3 0 .24 5
        ######## arsenate adsorption, model B, Table 2
        # 2Chalk_aOH-0.667 + CaAsO4- = Chalk_a2AsO4Ca-0.334 + 2OH-
        # log_k 3.59
        # -cd_music  0 -4 0 .26 5
        # 2Chalk_aOH-0.667 + CaHAsO4 = Chalk_a2HAsO4Ca+0.666 + 2OH-
        # log_k 1.27
        # -cd_music  0 -3 0 .24 5
    END
    """
end # cd_music_calcite

function diffuse_layer_oil(;oil_logk=[5.5,-4.7572,-3.8263,-3.4781, -8.93, -5.85, -5.85],
    oil_dh=[30,0.0,1.17152,-8.42239, 0.0, 0.0, 0.0])
    return """
        SURFACE_MASTER_SPECIES
            Oil_n      Oil_nH
            Oil_c      Oil_cOH
            Oil_w      Oil_wOH
        SURFACE_SPECIES
            Oil_nH = Oil_nH
            log_k   0.0
            Oil_cOH = Oil_cOH
            log_k   0.0
            Oil_wOH = Oil_wOH
            log_k   0.0
            # Oil surface (LLNL database for acetate)
            #oil 1
            Oil_nH + H+ = Oil_nH2+
            log_k   $(oil_logk[1])
            delta_h $(oil_dh[1])
            #oil 2
            Oil_cOH = Oil_cO- + H+
            log_k   $(oil_logk[2])
            delta_h $(oil_dh[2])
            #oil 3
            Oil_cOH + Ca+2 = Oil_cOCa+ + H+
            log_k   $(oil_logk[3])
            delta_h $(oil_dh[3])
            #oil 4
            Oil_cOH + Mg+2 = Oil_cOMg+ + H+
            log_k   $(oil_logk[4])
            delta_h $(oil_dh[4])
            #oil 5  
            Oil_wOH = Oil_wO- + H+
            log_k   $(oil_logk[5])
            delta_h $(oil_dh[5])
            #oil 6
            Oil_wOH + Ca+2 = Oil_wOCa+ + H+
            log_k   $(oil_logk[6])
            delta_h $(oil_dh[6])
            #oil 7
            Oil_wOH + Mg+2 = Oil_wOMg+ + H+
            log_k   $(oil_logk[7])
            delta_h $(oil_dh[7])
        END
    """
end

"""
from the work of Jin Song in Hirasaki's group
"""
function diffuse_layer_calcite(;log_k=[0.30, 1.74, 1.62, 0.14, 0.5, 2.23, 1.0, 0.09, -0.64], 
    dh=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
return """
    SURFACE_MASTER_SPECIES
        Chalk_a    Chalk_aOH-0.75
        Chalk_c    Chalk_cOH+0.75
    SURFACE_SPECIES
        Chalk_aOH-0.75 = Chalk_aOH-0.75
        log_k   0.0
        Chalk_cOH+0.75 = Chalk_cOH+0.75
        log_k   0.0
        # calcium site
        Chalk_aOH-0.75 + H+ = Chalk_aOH2+0.25
        log_k   $(log_k[1]) # 0.30
        delta_h      $(dh[1])
        Chalk_aOH-0.75 + Ca+2 = Chalk_aOHCa+1.25
        log_k   $(log_k[2]) # 1.74
        delta_h      $(dh[2])
        Chalk_aOH-0.75 + Mg+2 = Chalk_aOHMg+1.25
        log_k   $(log_k[3]) # 1.62
        delta_h      $(dh[3])
        Chalk_aOH-0.75 + Na+ = Chalk_aOHNa+0.25
        log_k   $(log_k[4]) # 0.14
        delta_h      $(dh[4])
        # carbonate site
        Chalk_cOH+0.75 + OH- = Chalk_cO-0.25 + H2O
        log_k   $(log_k[5]) # 0.5
        delta_h      $(dh[5])
        Chalk_cOH+0.75 + CO3-2 = Chalk_cOHCO3-1.25
        log_k   $(log_k[6]) # 2.23
        delta_h      $(dh[6])
        Chalk_cOH+0.75 + SO4-2 = Chalk_cOHSO4-1.25
        log_k   $(log_k[7]) # 1.0
        delta_h      $(dh[7])
        Chalk_cOH+0.75 + HCO3- = Chalk_cOH2CO3-0.25
        log_k   $(log_k[8]) # 0.09
        delta_h      $(dh[8])
        Chalk_cOH+0.75 + Cl- = Chalk_cOHCl-0.25
        log_k   $(log_k[9]) # -0.64
        delta_h      $(dh[9])
    END
"""
end

function calcite_surface()

end # calcite_surface