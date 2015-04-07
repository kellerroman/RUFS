subroutine wr(line,level)
implicit none
character(len=*), intent(in) :: line
integer         , intent(in) :: level
integer :: laenge
integer :: offset
character(len = 19) :: para_str
integer :: l,n1,n2

if (level <= 1) then
   laenge = 100
   offset = 1
else if (level == 2) then
   laenge = 80
   offset = 5
else
   laenge = 60
   offset = 10
end if
l = len_trim(line)

n1 = (laenge-l-2) / 2
n2 = laenge - n1 - l - 2

write(para_str,'("(",I2.2,"X,",I2.2,"A,X,A,X,",I2.2,"A)")') offset,n1,n2

!write(*,*) para_str

write(*,para_str) ("=",l = 1,n1),trim(line),("=",l=1,n2)

end subroutine wr