module m_mkfldh
contains
subroutine mkfldh(fieldname,ix,iy,layer,ia,ib,ic,length,fieldhead)
   character(len=*)  fieldname
   character(len=*)  fieldhead
   write(fieldhead,100) fieldname,layer,ix,iy,ia,ib,ic,length 
100 format("++ ",A5,3I2,3I4,I6)
end subroutine mkfldh
end module m_mkfldh
