module commonstuff_hexa

integer*8 locedg_hexa(12,2)
data locedg_hexa /1,2,3,4,5,6,7,8,1,2,3,4, &
             2,3,4,1,6,7,8,5,5,6,7,8/

integer*8 locfacedg_hexa(4,2)
data locfacedg_hexa /2,3,4,1, &
                3,4,1,2/

integer*8 locfac_hexa(6,4)
data locfac_hexa /1,4,8,5,2,4, &
             2,3,7,6,6,8, &
             3,7,6,2,7,5, &	 
			 4,8,5,1,3,1/

logical robust

end module
