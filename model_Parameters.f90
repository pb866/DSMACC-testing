! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Parameter Module File
! 
! Generated by KPP-2.2.4_gc symbolic chemistry Kinetics PreProcessor
!       (http://www.cs.vt.edu/~asandu/Software/KPP)
! KPP is distributed under GPL, the general public licence
!       (http://www.gnu.org/copyleft/gpl.html)
! (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa
! (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech
!     With important contributions from:
!        M. Damian, Villanova University, USA
!        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
! 
! File                 : model_Parameters.f90
! Time                 : Sun Feb  5 23:53:52 2017
! Working directory    : /work/home/dp626/DSMACC-testing
! Equation file        : model.kpp
! Output root filename : model
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE model_Parameters

  USE model_Precision
  PUBLIC
  SAVE


! NSPEC - Number of chemical species
  INTEGER, PARAMETER :: NSPEC = 456 
! NVAR - Number of Variable species
  INTEGER, PARAMETER :: NVAR = 455 
! NFLUX - Number of Reaction Flux species
  INTEGER, PARAMETER :: NFLUX = 1 
! NFAM - Number of Prod/Loss Families
  INTEGER, PARAMETER :: NFAM = 1 
! NVARACT - Number of Active species
  INTEGER, PARAMETER :: NVARACT = 454 
! NFIX - Number of Fixed species
  INTEGER, PARAMETER :: NFIX = 1 
! NREACT - Number of reactions
  INTEGER, PARAMETER :: NREACT = 1452 
! NVARST - Starting of variables in conc. vect.
  INTEGER, PARAMETER :: NVARST = 1 
! NFIXST - Starting of fixed in conc. vect.
  INTEGER, PARAMETER :: NFIXST = 456 
! NONZERO - Number of nonzero entries in Jacobian
  INTEGER, PARAMETER :: NONZERO = 4118 
! LU_NONZERO - Number of nonzero entries in LU factoriz. of Jacobian
  INTEGER, PARAMETER :: LU_NONZERO = 5137 
! CNVAR - (NVAR+1) Number of elements in compressed row format
  INTEGER, PARAMETER :: CNVAR = 456 
! NLOOKAT - Number of species to look at
  INTEGER, PARAMETER :: NLOOKAT = 0 
! NMONITOR - Number of species to monitor
  INTEGER, PARAMETER :: NMONITOR = 0 
! NMASS - Number of atoms to check mass balance
  INTEGER, PARAMETER :: NMASS = 1 

! Index declaration for variable species in C and VAR
!   VAR(ind_spc) = C(ind_spc)

  INTEGER, PARAMETER :: ind_DUMMY = 1 
  INTEGER, PARAMETER :: ind_SA = 2 
  INTEGER, PARAMETER :: ind_SO3 = 3 
  INTEGER, PARAMETER :: ind_CL = 4 
  INTEGER, PARAMETER :: ind_O1D = 5 
  INTEGER, PARAMETER :: ind_HSO3 = 6 
  INTEGER, PARAMETER :: ind_PROPALO = 7 
  INTEGER, PARAMETER :: ind_INB1HPCO2H = 8 
  INTEGER, PARAMETER :: ind_BIACETO = 9 
  INTEGER, PARAMETER :: ind_C510OH = 10 
  INTEGER, PARAMETER :: ind_INANCOCHO = 11 
  INTEGER, PARAMETER :: ind_CH4 = 12 
  INTEGER, PARAMETER :: ind_HCOOH = 13 
  INTEGER, PARAMETER :: ind_CH3O2NO2 = 14 
  INTEGER, PARAMETER :: ind_C57OH = 15 
  INTEGER, PARAMETER :: ind_ISOPBOH = 16 
  INTEGER, PARAMETER :: ind_C58OH = 17 
  INTEGER, PARAMETER :: ind_N2O5 = 18 
  INTEGER, PARAMETER :: ind_IEPOXA = 19 
  INTEGER, PARAMETER :: ind_ETHGLY = 20 
  INTEGER, PARAMETER :: ind_INAOH = 21 
  INTEGER, PARAMETER :: ind_INCOH = 22 
  INTEGER, PARAMETER :: ind_HIEPOXB = 23 
  INTEGER, PARAMETER :: ind_ALLYLOH = 24 
  INTEGER, PARAMETER :: ind_INDOH = 25 
  INTEGER, PARAMETER :: ind_C524OH = 26 
  INTEGER, PARAMETER :: ind_OCCOHCOH = 27 
  INTEGER, PARAMETER :: ind_C2OHOCO2H = 28 
  INTEGER, PARAMETER :: ind_ISOPDOH = 29 
  INTEGER, PARAMETER :: ind_CH3OH = 30 
  INTEGER, PARAMETER :: ind_H2 = 31 
  INTEGER, PARAMETER :: ind_INB2O = 32 
  INTEGER, PARAMETER :: ind_HOCH2COCO2H = 33 
  INTEGER, PARAMETER :: ind_NC524OH = 34 
  INTEGER, PARAMETER :: ind_HONO = 35 
  INTEGER, PARAMETER :: ind_C59O = 36 
  INTEGER, PARAMETER :: ind_IEPOXC = 37 
  INTEGER, PARAMETER :: ind_MACO2H = 38 
  INTEGER, PARAMETER :: ind_MGLYOOB = 39 
  INTEGER, PARAMETER :: ind_MGLYOOA = 40 
  INTEGER, PARAMETER :: ind_H1CO23CHO = 41 
  INTEGER, PARAMETER :: ind_MVKO = 42 
  INTEGER, PARAMETER :: ind_CH2OOG = 43 
  INTEGER, PARAMETER :: ind_CH2OOE = 44 
  INTEGER, PARAMETER :: ind_CH2OOC = 45 
  INTEGER, PARAMETER :: ind_NC4OOA = 46 
  INTEGER, PARAMETER :: ind_NOAOOA = 47 
  INTEGER, PARAMETER :: ind_C525O = 48 
  INTEGER, PARAMETER :: ind_HPNC524CO = 49 
  INTEGER, PARAMETER :: ind_HMVKBO = 50 
  INTEGER, PARAMETER :: ind_HMGLYOOA = 51 
  INTEGER, PARAMETER :: ind_HNC524CO = 52 
  INTEGER, PARAMETER :: ind_A2PANO = 53 
  INTEGER, PARAMETER :: ind_MVKOHBO = 54 
  INTEGER, PARAMETER :: ind_NC3OOA = 55 
  INTEGER, PARAMETER :: ind_CH3CO2H = 56 
  INTEGER, PARAMETER :: ind_NC2OOA = 57 
  INTEGER, PARAMETER :: ind_INAHPCO3H = 58 
  INTEGER, PARAMETER :: ind_HC4CCO2H = 59 
  INTEGER, PARAMETER :: ind_CH3CO3H = 60 
  INTEGER, PARAMETER :: ind_IECCO3H = 61 
  INTEGER, PARAMETER :: ind_NC524NO3 = 62 
  INTEGER, PARAMETER :: ind_OCCOHCOOH = 63 
  INTEGER, PARAMETER :: ind_C58NO3CO2H = 64 
  INTEGER, PARAMETER :: ind_INB1NACO2H = 65 
  INTEGER, PARAMETER :: ind_MVKOHAOOH = 66 
  INTEGER, PARAMETER :: ind_HOCHOCOOH = 67 
  INTEGER, PARAMETER :: ind_IEPOXB = 68 
  INTEGER, PARAMETER :: ind_C2OHOCOOH = 69 
  INTEGER, PARAMETER :: ind_INB1HPCO3H = 70 
  INTEGER, PARAMETER :: ind_C510OOH = 71 
  INTEGER, PARAMETER :: ind_INB1NBCO2H = 72 
  INTEGER, PARAMETER :: ind_MACO3H = 73 
  INTEGER, PARAMETER :: ind_INCNCO2H = 74 
  INTEGER, PARAMETER :: ind_C57OOH = 75 
  INTEGER, PARAMETER :: ind_NISOPOOH = 76 
  INTEGER, PARAMETER :: ind_IPRHOCO2H = 77 
  INTEGER, PARAMETER :: ind_MVKOHANO3 = 78 
  INTEGER, PARAMETER :: ind_HMACO2H = 79 
  INTEGER, PARAMETER :: ind_IEACO3H = 80 
  INTEGER, PARAMETER :: ind_C3DIOLOOH = 81 
  INTEGER, PARAMETER :: ind_C57NO3CO2H = 82 
  INTEGER, PARAMETER :: ind_HMVKAOOH = 83 
  INTEGER, PARAMETER :: ind_INAHPCO2H = 84 
  INTEGER, PARAMETER :: ind_INAHCO2H = 85 
  INTEGER, PARAMETER :: ind_PRNO3CO2H = 86 
  INTEGER, PARAMETER :: ind_INANCO2H = 87 
  INTEGER, PARAMETER :: ind_INDHPCO3H = 88 
  INTEGER, PARAMETER :: ind_C58OOH = 89 
  INTEGER, PARAMETER :: ind_HO2NO2 = 90 
  INTEGER, PARAMETER :: ind_INAHCHO = 91 
  INTEGER, PARAMETER :: ind_IBUTALOH = 92 
  INTEGER, PARAMETER :: ind_INB1HPCHO = 93 
  INTEGER, PARAMETER :: ind_CH3CHOOA = 94 
  INTEGER, PARAMETER :: ind_CH3O = 95 
  INTEGER, PARAMETER :: ind_HMVKNGLYOX = 96 
  INTEGER, PARAMETER :: ind_INDHCHO = 97 
  INTEGER, PARAMETER :: ind_INDHPCHO = 98 
  INTEGER, PARAMETER :: ind_INAHPCHO = 99 
  INTEGER, PARAMETER :: ind_GLYOOB = 100 
  INTEGER, PARAMETER :: ind_INCCO = 101 
  INTEGER, PARAMETER :: ind_INB1O = 102 
  INTEGER, PARAMETER :: ind_CO2N3CHO = 103 
  INTEGER, PARAMETER :: ind_HYPROPO2H = 104 
  INTEGER, PARAMETER :: ind_ACO2H = 105 
  INTEGER, PARAMETER :: ind_INB1GLYOX = 106 
  INTEGER, PARAMETER :: ind_C524CO = 107 
  INTEGER, PARAMETER :: ind_ETHO2HNO3 = 108 
  INTEGER, PARAMETER :: ind_IEB1O = 109 
  INTEGER, PARAMETER :: ind_IPROPOLPER = 110 
  INTEGER, PARAMETER :: ind_IPROPOLO = 111 
  INTEGER, PARAMETER :: ind_CH2OOA = 112 
  INTEGER, PARAMETER :: ind_HOCH2CO2H = 113 
  INTEGER, PARAMETER :: ind_IPRHOCO3H = 114 
  INTEGER, PARAMETER :: ind_INB1NBCHO = 115 
  INTEGER, PARAMETER :: ind_CH3OOH = 116 
  INTEGER, PARAMETER :: ind_PROPGLY = 117 
  INTEGER, PARAMETER :: ind_INCO = 118 
  INTEGER, PARAMETER :: ind_NISOPNO3 = 119 
  INTEGER, PARAMETER :: ind_NISOPO = 120 
  INTEGER, PARAMETER :: ind_IEC2OOH = 121 
  INTEGER, PARAMETER :: ind_H13CO2CO3H = 122 
  INTEGER, PARAMETER :: ind_HMGLOOA = 123 
  INTEGER, PARAMETER :: ind_CH3COCO2H = 124 
  INTEGER, PARAMETER :: ind_INB1NACO3H = 125 
  INTEGER, PARAMETER :: ind_HIEB1O = 126 
  INTEGER, PARAMETER :: ind_INB1HPPAN = 127 
  INTEGER, PARAMETER :: ind_CH3NO3 = 128 
  INTEGER, PARAMETER :: ind_ISOPAOH = 129 
  INTEGER, PARAMETER :: ind_INB1NBCO3H = 130 
  INTEGER, PARAMETER :: ind_HYETHO2H = 131 
  INTEGER, PARAMETER :: ind_ETHOHNO3 = 132 
  INTEGER, PARAMETER :: ind_H14CO23C4 = 133 
  INTEGER, PARAMETER :: ind_CO2N3CO3H = 134 
  INTEGER, PARAMETER :: ind_INDHCO3H = 135 
  INTEGER, PARAMETER :: ind_INB1NACHO = 136 
  INTEGER, PARAMETER :: ind_MACRNBPAN = 137 
  INTEGER, PARAMETER :: ind_IEB2O = 138 
  INTEGER, PARAMETER :: ind_HC4ACO2H = 139 
  INTEGER, PARAMETER :: ind_NO3CH2CO3H = 140 
  INTEGER, PARAMETER :: ind_HMACO3H = 141 
  INTEGER, PARAMETER :: ind_IEC1O = 142 
  INTEGER, PARAMETER :: ind_HOCH2CO3H = 143 
  INTEGER, PARAMETER :: ind_PR2O2HNO3 = 144 
  INTEGER, PARAMETER :: ind_INAHCO3H = 145 
  INTEGER, PARAMETER :: ind_HCOCOHCO3H = 146 
  INTEGER, PARAMETER :: ind_CH3COCH2O = 147 
  INTEGER, PARAMETER :: ind_ACLOOA = 148 
  INTEGER, PARAMETER :: ind_INAO = 149 
  INTEGER, PARAMETER :: ind_INANCHO = 150 
  INTEGER, PARAMETER :: ind_C524NO3 = 151 
  INTEGER, PARAMETER :: ind_MACRO = 152 
  INTEGER, PARAMETER :: ind_MACROOH = 153 
  INTEGER, PARAMETER :: ind_GLYOOC = 154 
  INTEGER, PARAMETER :: ind_HMVKBOOH = 155 
  INTEGER, PARAMETER :: ind_HIEB2O = 156 
  INTEGER, PARAMETER :: ind_PR1O2HNO3 = 157 
  INTEGER, PARAMETER :: ind_C59OOH = 158 
  INTEGER, PARAMETER :: ind_INCNCO3H = 159 
  INTEGER, PARAMETER :: ind_PROLNO3 = 160 
  INTEGER, PARAMETER :: ind_HC4CCO3H = 161 
  INTEGER, PARAMETER :: ind_OCCOHCO = 162 
  INTEGER, PARAMETER :: ind_PRNO3CO3H = 163 
  INTEGER, PARAMETER :: ind_NO3CH2CO2H = 164 
  INTEGER, PARAMETER :: ind_INB1OH = 165 
  INTEGER, PARAMETER :: ind_CHOMOHCO3H = 166 
  INTEGER, PARAMETER :: ind_NC4CO2H = 167 
  INTEGER, PARAMETER :: ind_INANCO3H = 168 
  INTEGER, PARAMETER :: ind_C525OOH = 169 
  INTEGER, PARAMETER :: ind_MGLOOA = 170 
  INTEGER, PARAMETER :: ind_MACRNOOA = 171 
  INTEGER, PARAMETER :: ind_INCNCHO = 172 
  INTEGER, PARAMETER :: ind_MACROOA = 173 
  INTEGER, PARAMETER :: ind_C510O = 174 
  INTEGER, PARAMETER :: ind_C57O = 175 
  INTEGER, PARAMETER :: ind_C524O = 176 
  INTEGER, PARAMETER :: ind_PHAN = 177 
  INTEGER, PARAMETER :: ind_IEAPAN = 178 
  INTEGER, PARAMETER :: ind_HMVKAO = 179 
  INTEGER, PARAMETER :: ind_PAN = 180 
  INTEGER, PARAMETER :: ind_ISOPDO = 181 
  INTEGER, PARAMETER :: ind_DNC524CO = 182 
  INTEGER, PARAMETER :: ind_CO2N3PAN = 183 
  INTEGER, PARAMETER :: ind_IECPAN = 184 
  INTEGER, PARAMETER :: ind_ACRPAN = 185 
  INTEGER, PARAMETER :: ind_HMACROH = 186 
  INTEGER, PARAMETER :: ind_ISOPAO = 187 
  INTEGER, PARAMETER :: ind_GAOOB = 188 
  INTEGER, PARAMETER :: ind_HYPROPO = 189 
  INTEGER, PARAMETER :: ind_ISOPBOOH = 190 
  INTEGER, PARAMETER :: ind_IPROPOLPAN = 191 
  INTEGER, PARAMETER :: ind_C5PAN19 = 192 
  INTEGER, PARAMETER :: ind_C5PAN18 = 193 
  INTEGER, PARAMETER :: ind_C5PAN17 = 194 
  INTEGER, PARAMETER :: ind_HCOC5 = 195 
  INTEGER, PARAMETER :: ind_PRNO3PAN = 196 
  INTEGER, PARAMETER :: ind_HCOCH2O = 197 
  INTEGER, PARAMETER :: ind_ISOPCO = 198 
  INTEGER, PARAMETER :: ind_HOCH2CH2O = 199 
  INTEGER, PARAMETER :: ind_INDHPPAN = 200 
  INTEGER, PARAMETER :: ind_INB2OOH = 201 
  INTEGER, PARAMETER :: ind_C58NO3PAN = 202 
  INTEGER, PARAMETER :: ind_NC4CO3H = 203 
  INTEGER, PARAMETER :: ind_C58NO3CO3H = 204 
  INTEGER, PARAMETER :: ind_INB1NAPAN = 205 
  INTEGER, PARAMETER :: ind_INCNPAN = 206 
  INTEGER, PARAMETER :: ind_HCOCOHPAN = 207 
  INTEGER, PARAMETER :: ind_CH3COCO3H = 208 
  INTEGER, PARAMETER :: ind_C32OH13CO = 209 
  INTEGER, PARAMETER :: ind_INANPAN = 210 
  INTEGER, PARAMETER :: ind_MACROHO = 211 
  INTEGER, PARAMETER :: ind_C57NO3PAN = 212 
  INTEGER, PARAMETER :: ind_C58O = 213 
  INTEGER, PARAMETER :: ind_CHOMOHPAN = 214 
  INTEGER, PARAMETER :: ind_INDHPAN = 215 
  INTEGER, PARAMETER :: ind_INB1NBPAN = 216 
  INTEGER, PARAMETER :: ind_HC4ACO3H = 217 
  INTEGER, PARAMETER :: ind_HYPERACET = 218 
  INTEGER, PARAMETER :: ind_MVKOOH = 219 
  INTEGER, PARAMETER :: ind_MVKOHAO = 220 
  INTEGER, PARAMETER :: ind_A2PAN = 221 
  INTEGER, PARAMETER :: ind_MACRNBCO2H = 222 
  INTEGER, PARAMETER :: ind_MVKOHBOOH = 223 
  INTEGER, PARAMETER :: ind_C4PAN10 = 224 
  INTEGER, PARAMETER :: ind_NO3CH2PAN = 225 
  INTEGER, PARAMETER :: ind_MACRNPAN = 226 
  INTEGER, PARAMETER :: ind_ACO3H = 227 
  INTEGER, PARAMETER :: ind_C3DIOLO = 228 
  INTEGER, PARAMETER :: ind_CHOCOHCO = 229 
  INTEGER, PARAMETER :: ind_C4PAN5 = 230 
  INTEGER, PARAMETER :: ind_C4PAN6 = 231 
  INTEGER, PARAMETER :: ind_INANCO = 232 
  INTEGER, PARAMETER :: ind_C4M2ALOHNO3 = 233 
  INTEGER, PARAMETER :: ind_MVKOOA = 234 
  INTEGER, PARAMETER :: ind_INCGLYOX = 235 
  INTEGER, PARAMETER :: ind_CO2H3CO3H = 236 
  INTEGER, PARAMETER :: ind_HIEB1OOH = 237 
  INTEGER, PARAMETER :: ind_IEB2OOH = 238 
  INTEGER, PARAMETER :: ind_PROPOLNO3 = 239 
  INTEGER, PARAMETER :: ind_INAHPPAN = 240 
  INTEGER, PARAMETER :: ind_IEB1OOH = 241 
  INTEGER, PARAMETER :: ind_HIEB2OOH = 242 
  INTEGER, PARAMETER :: ind_HCOCH2OOH = 243 
  INTEGER, PARAMETER :: ind_NC524OOH = 244 
  INTEGER, PARAMETER :: ind_INANCOCO2H = 245 
  INTEGER, PARAMETER :: ind_HCOCO2H = 246 
  INTEGER, PARAMETER :: ind_COHM2CO2H = 247 
  INTEGER, PARAMETER :: ind_IPROPOLO2H = 248 
  INTEGER, PARAMETER :: ind_MACRNBCO3H = 249 
  INTEGER, PARAMETER :: ind_HMACROOH = 250 
  INTEGER, PARAMETER :: ind_HMACRO = 251 
  INTEGER, PARAMETER :: ind_C57NO3CO3H = 252 
  INTEGER, PARAMETER :: ind_INAHPAN = 253 
  INTEGER, PARAMETER :: ind_MACRNCO3H = 254 
  INTEGER, PARAMETER :: ind_HMPAN = 255 
  INTEGER, PARAMETER :: ind_ISOPBO = 256 
  INTEGER, PARAMETER :: ind_C42AOH = 257 
  INTEGER, PARAMETER :: ind_ETHENO3O = 258 
  INTEGER, PARAMETER :: ind_GLYPAN = 259 
  INTEGER, PARAMETER :: ind_MACRNCO2H = 260 
  INTEGER, PARAMETER :: ind_C3MDIALOH = 261 
  INTEGER, PARAMETER :: ind_INCOOH = 262 
  INTEGER, PARAMETER :: ind_INANCOPAN = 263 
  INTEGER, PARAMETER :: ind_ISOPAOOH = 264 
  INTEGER, PARAMETER :: ind_IECCHO = 265 
  INTEGER, PARAMETER :: ind_MACROHNO3 = 266 
  INTEGER, PARAMETER :: ind_INDOOH = 267 
  INTEGER, PARAMETER :: ind_INANCOCO3H = 268 
  INTEGER, PARAMETER :: ind_MACROHOOH = 269 
  INTEGER, PARAMETER :: ind_COHM2CO3H = 270 
  INTEGER, PARAMETER :: ind_INAOOH = 271 
  INTEGER, PARAMETER :: ind_MMALNACO2H = 272 
  INTEGER, PARAMETER :: ind_MMALNBCO2H = 273 
  INTEGER, PARAMETER :: ind_VGLYOX = 274 
  INTEGER, PARAMETER :: ind_CH2OOB = 275 
  INTEGER, PARAMETER :: ind_PRONO3BO = 276 
  INTEGER, PARAMETER :: ind_HMVKNO3 = 277 
  INTEGER, PARAMETER :: ind_CH3COPAN = 278 
  INTEGER, PARAMETER :: ind_O = 279 
  INTEGER, PARAMETER :: ind_MVKOHAOH = 280 
  INTEGER, PARAMETER :: ind_PRONO3AO = 281 
  INTEGER, PARAMETER :: ind_IEC1OOH = 282 
  INTEGER, PARAMETER :: ind_INCNO3 = 283 
  INTEGER, PARAMETER :: ind_C524OOH = 284 
  INTEGER, PARAMETER :: ind_MMALNACO3H = 285 
  INTEGER, PARAMETER :: ind_MMALNBCO3H = 286 
  INTEGER, PARAMETER :: ind_IEACHO = 287 
  INTEGER, PARAMETER :: ind_CO23C3CHO = 288 
  INTEGER, PARAMETER :: ind_CH3COCH3 = 289 
  INTEGER, PARAMETER :: ind_NISOPO2 = 290 
  INTEGER, PARAMETER :: ind_INB1CO = 291 
  INTEGER, PARAMETER :: ind_CH3C2H2O2 = 292 
  INTEGER, PARAMETER :: ind_ISOPDOOH = 293 
  INTEGER, PARAMETER :: ind_INDO = 294 
  INTEGER, PARAMETER :: ind_MMALNBPAN = 295 
  INTEGER, PARAMETER :: ind_INB1NO3 = 296 
  INTEGER, PARAMETER :: ind_NMGLYOX = 297 
  INTEGER, PARAMETER :: ind_HCOCO3H = 298 
  INTEGER, PARAMETER :: ind_COHM2PAN = 299 
  INTEGER, PARAMETER :: ind_HMVKANO3 = 300 
  INTEGER, PARAMETER :: ind_INB1OOH = 301 
  INTEGER, PARAMETER :: ind_ISOPCOOH = 302 
  INTEGER, PARAMETER :: ind_C3DIOLO2 = 303 
  INTEGER, PARAMETER :: ind_C58NO3 = 304 
  INTEGER, PARAMETER :: ind_HMACRO2 = 305 
  INTEGER, PARAMETER :: ind_NC524O = 306 
  INTEGER, PARAMETER :: ind_HOCH2CH2O2 = 307 
  INTEGER, PARAMETER :: ind_C525O2 = 308 
  INTEGER, PARAMETER :: ind_C510O2 = 309 
  INTEGER, PARAMETER :: ind_C57O2 = 310 
  INTEGER, PARAMETER :: ind_MACRNB = 311 
  INTEGER, PARAMETER :: ind_C59O2 = 312 
  INTEGER, PARAMETER :: ind_MPAN = 313 
  INTEGER, PARAMETER :: ind_INDHPCO3 = 314 
  INTEGER, PARAMETER :: ind_INDHCO3 = 315 
  INTEGER, PARAMETER :: ind_CO2N3CO3 = 316 
  INTEGER, PARAMETER :: ind_H2O2 = 317 
  INTEGER, PARAMETER :: ind_MMALNAPAN = 318 
  INTEGER, PARAMETER :: ind_MVKOH = 319 
  INTEGER, PARAMETER :: ind_ETHENO3O2 = 320 
  INTEGER, PARAMETER :: ind_C57NO3 = 321 
  INTEGER, PARAMETER :: ind_HIEB2O2 = 322 
  INTEGER, PARAMETER :: ind_CHOPRNO3 = 323 
  INTEGER, PARAMETER :: ind_MACROH = 324 
  INTEGER, PARAMETER :: ind_CONM2CO2H = 325 
  INTEGER, PARAMETER :: ind_H13CO2CO3 = 326 
  INTEGER, PARAMETER :: ind_IEACO3 = 327 
  INTEGER, PARAMETER :: ind_HCOCOHCO3 = 328 
  INTEGER, PARAMETER :: ind_IECCO3 = 329 
  INTEGER, PARAMETER :: ind_CH3CHOHCO3 = 330 
  INTEGER, PARAMETER :: ind_C2H4 = 331 
  INTEGER, PARAMETER :: ind_PRONO3BO2 = 332 
  INTEGER, PARAMETER :: ind_HMVKAO2 = 333 
  INTEGER, PARAMETER :: ind_HYPROPO2 = 334 
  INTEGER, PARAMETER :: ind_H13CO2C3 = 335 
  INTEGER, PARAMETER :: ind_ACO3B = 336 
  INTEGER, PARAMETER :: ind_IPROPOLO2 = 337 
  INTEGER, PARAMETER :: ind_CH3CHOHCHO = 338 
  INTEGER, PARAMETER :: ind_HCOCH2O2 = 339 
  INTEGER, PARAMETER :: ind_HMGLYOO = 340 
  INTEGER, PARAMETER :: ind_NC4OO = 341 
  INTEGER, PARAMETER :: ind_MVKOO = 342 
  INTEGER, PARAMETER :: ind_CONM2CO3H = 343 
  INTEGER, PARAMETER :: ind_NC4CO3 = 344 
  INTEGER, PARAMETER :: ind_INB1HPCO3 = 345 
  INTEGER, PARAMETER :: ind_NC3OO = 346 
  INTEGER, PARAMETER :: ind_PRNO3CO3 = 347 
  INTEGER, PARAMETER :: ind_NOAOO = 348 
  INTEGER, PARAMETER :: ind_INANCO3 = 349 
  INTEGER, PARAMETER :: ind_INCO2 = 350 
  INTEGER, PARAMETER :: ind_ISOPBO2 = 351 
  INTEGER, PARAMETER :: ind_C58O2 = 352 
  INTEGER, PARAMETER :: ind_NC524O2 = 353 
  INTEGER, PARAMETER :: ind_HIEB1O2 = 354 
  INTEGER, PARAMETER :: ind_ISOPANO3 = 355 
  INTEGER, PARAMETER :: ind_INAHCO3 = 356 
  INTEGER, PARAMETER :: ind_INAHPCO3 = 357 
  INTEGER, PARAMETER :: ind_CH3CHOO = 358 
  INTEGER, PARAMETER :: ind_MVKOHBO2 = 359 
  INTEGER, PARAMETER :: ind_CONM2PAN = 360 
  INTEGER, PARAMETER :: ind_CO2H3CO3 = 361 
  INTEGER, PARAMETER :: ind_IEB1O2 = 362 
  INTEGER, PARAMETER :: ind_INCNCO3 = 363 
  INTEGER, PARAMETER :: ind_CONM2CHO = 364 
  INTEGER, PARAMETER :: ind_MGLOO = 365 
  INTEGER, PARAMETER :: ind_COHM2CO3 = 366 
  INTEGER, PARAMETER :: ind_PRONO3AO2 = 367 
  INTEGER, PARAMETER :: ind_IEC1O2 = 368 
  INTEGER, PARAMETER :: ind_INB1NACO3 = 369 
  INTEGER, PARAMETER :: ind_INB1NBCO3 = 370 
  INTEGER, PARAMETER :: ind_MGLYOO = 371 
  INTEGER, PARAMETER :: ind_ISOPDNO3 = 372 
  INTEGER, PARAMETER :: ind_MVKOHAO2 = 373 
  INTEGER, PARAMETER :: ind_CO23C4NO3 = 374 
  INTEGER, PARAMETER :: ind_MACRNBCO3 = 375 
  INTEGER, PARAMETER :: ind_HC4CCO3 = 376 
  INTEGER, PARAMETER :: ind_NOA = 377 
  INTEGER, PARAMETER :: ind_IEB2O2 = 378 
  INTEGER, PARAMETER :: ind_HMACO3 = 379 
  INTEGER, PARAMETER :: ind_HC4ACO3 = 380 
  INTEGER, PARAMETER :: ind_IPRHOCO3 = 381 
  INTEGER, PARAMETER :: ind_CH3COCH2O2 = 382 
  INTEGER, PARAMETER :: ind_HO12CO3C4 = 383 
  INTEGER, PARAMETER :: ind_HMVKBO2 = 384 
  INTEGER, PARAMETER :: ind_A2PANOO = 385 
  INTEGER, PARAMETER :: ind_NC2OO = 386 
  INTEGER, PARAMETER :: ind_GAOO = 387 
  INTEGER, PARAMETER :: ind_HMGLOO = 388 
  INTEGER, PARAMETER :: ind_H13CO2CHO = 389 
  INTEGER, PARAMETER :: ind_ACO3 = 390 
  INTEGER, PARAMETER :: ind_MACRO2 = 391 
  INTEGER, PARAMETER :: ind_ACR = 392 
  INTEGER, PARAMETER :: ind_ACLOO = 393 
  INTEGER, PARAMETER :: ind_INANO3 = 394 
  INTEGER, PARAMETER :: ind_MACROHO2 = 395 
  INTEGER, PARAMETER :: ind_HOCH2COCHO = 396 
  INTEGER, PARAMETER :: ind_MACRNOO = 397 
  INTEGER, PARAMETER :: ind_MACROO = 398 
  INTEGER, PARAMETER :: ind_MMALNACO3 = 399 
  INTEGER, PARAMETER :: ind_CHOMOHCO3 = 400 
  INTEGER, PARAMETER :: ind_INB1O2 = 401 
  INTEGER, PARAMETER :: ind_C524O2 = 402 
  INTEGER, PARAMETER :: ind_INDO2 = 403 
  INTEGER, PARAMETER :: ind_HMACR = 404 
  INTEGER, PARAMETER :: ind_ISOPCO2 = 405 
  INTEGER, PARAMETER :: ind_MACO3 = 406 
  INTEGER, PARAMETER :: ind_CH3O2 = 407 
  INTEGER, PARAMETER :: ind_C3H6 = 408 
  INTEGER, PARAMETER :: ind_HNO3 = 409 
  INTEGER, PARAMETER :: ind_OCCOHCO2 = 410 
  INTEGER, PARAMETER :: ind_NC4CHO = 411 
  INTEGER, PARAMETER :: ind_CH3COCO3 = 412 
  INTEGER, PARAMETER :: ind_INANCOCO3 = 413 
  INTEGER, PARAMETER :: ind_NO3CH2CO3 = 414 
  INTEGER, PARAMETER :: ind_HCOCO3 = 415 
  INTEGER, PARAMETER :: ind_MMALNBCO3 = 416 
  INTEGER, PARAMETER :: ind_CO2H3CHO = 417 
  INTEGER, PARAMETER :: ind_INB2O2 = 418 
  INTEGER, PARAMETER :: ind_C58NO3CO3 = 419 
  INTEGER, PARAMETER :: ind_C57NO3CO3 = 420 
  INTEGER, PARAMETER :: ind_C5H8 = 421 
  INTEGER, PARAMETER :: ind_ISOPCNO3 = 422 
  INTEGER, PARAMETER :: ind_MVKO2 = 423 
  INTEGER, PARAMETER :: ind_GLYOO = 424 
  INTEGER, PARAMETER :: ind_HC4CCHO = 425 
  INTEGER, PARAMETER :: ind_HC4ACHO = 426 
  INTEGER, PARAMETER :: ind_ACETOL = 427 
  INTEGER, PARAMETER :: ind_INAO2 = 428 
  INTEGER, PARAMETER :: ind_NO3CH2CHO = 429 
  INTEGER, PARAMETER :: ind_CH3CHO = 430 
  INTEGER, PARAMETER :: ind_ISOPDO2 = 431 
  INTEGER, PARAMETER :: ind_GLYOX = 432 
  INTEGER, PARAMETER :: ind_MVKNO3 = 433 
  INTEGER, PARAMETER :: ind_MGLYOX = 434 
  INTEGER, PARAMETER :: ind_ISOPAO2 = 435 
  INTEGER, PARAMETER :: ind_MACRNCO3 = 436 
  INTEGER, PARAMETER :: ind_CONM2CO3 = 437 
  INTEGER, PARAMETER :: ind_ISOPBNO3 = 438 
  INTEGER, PARAMETER :: ind_MACRNO3 = 439 
  INTEGER, PARAMETER :: ind_HCHO = 440 
  INTEGER, PARAMETER :: ind_CH2OO = 441 
  INTEGER, PARAMETER :: ind_HOCH2CO3 = 442 
  INTEGER, PARAMETER :: ind_HOCH2CHO = 443 
  INTEGER, PARAMETER :: ind_MACR = 444 
  INTEGER, PARAMETER :: ind_CH3CO3 = 445 
  INTEGER, PARAMETER :: ind_BIACETOH = 446 
  INTEGER, PARAMETER :: ind_HO2 = 447 
  INTEGER, PARAMETER :: ind_SO2 = 448 
  INTEGER, PARAMETER :: ind_NO = 449 
  INTEGER, PARAMETER :: ind_CO = 450 
  INTEGER, PARAMETER :: ind_OH = 451 
  INTEGER, PARAMETER :: ind_MVK = 452 
  INTEGER, PARAMETER :: ind_O3 = 453 
  INTEGER, PARAMETER :: ind_NO2 = 454 
  INTEGER, PARAMETER :: ind_NO3 = 455 

! Index declaration for fixed species in C
!   C(ind_spc)

  INTEGER, PARAMETER :: ind_EMISS = 456 

! Index declaration for dummy species

  INTEGER, PARAMETER :: ind_H = 0 
  INTEGER, PARAMETER :: ind_O2 = 0 
  INTEGER, PARAMETER :: ind_A = 0 
  INTEGER, PARAMETER :: ind_O3FT = 0 
  INTEGER, PARAMETER :: ind_NA = 0 

! Index declaration for fixed species in FIX
!    FIX(indf_spc) = C(ind_spc) = C(NVAR+indf_spc)

  INTEGER, PARAMETER :: indf_EMISS = 1 

END MODULE model_Parameters

