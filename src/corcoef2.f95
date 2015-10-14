subroutine corcoef(tmpm,N,M,cormatrix)

!real, dimension(N,M):: tmpm
!real, dimension(M,M) :: cor_matrix
!real:: output=0

integer N,M
double precision tmpm(N,M), cormatrix(M,M), output
output = 0.0

do i = 1,M
        do j = 1,M
                if (i == j) then
                        output = 1.0
                        cormatrix(i,j) = output
                else
                        call cor(N, tmpm(:,i), tmpm(:,j), output)
                        cormatrix(i,j) = output
                end if
        end do
end do

return

end subroutine corcoef
