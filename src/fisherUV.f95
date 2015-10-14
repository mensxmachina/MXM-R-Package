subroutine fisher_uv(R, C, y, dataset, cs_cols, pvalues, stats, targetID)

        integer R,C,cs_cols,e,targetID
        integer, dimension(2):: xyIdx
        integer, dimension(cs_cols):: csIdx
        double precision dataset(R,C)
        integer xIndex, flag
        double precision y(R)
        integer csIndex(cs_cols)
        double precision pvalue, stat, fn_val
        double precision invm(cs_cols,cs_cols)
        real, dimension(2,2):: rcm
        double precision cor_matrix(cs_cols+2,cs_cols+2)

        double precision tmpm(R,cs_cols+2)
        real, dimension(size(dataset,1),size(csIndex)):: cs
        double precision x(R)
		double precision df, w
		
        real, dimension(1,2):: a
        real:: z
        real:: t1,t2,t3

		double complex w_complex;
        complex z_complex;
		
		double precision pvalues(C), stats(C)
		integer iter;

		cs_cols = 0;
		
		do iter=1,C
		
			xIndex = iter
      
      if( targetID == iter ) then
        pvalues(iter) = 1
        stats(iter) = 0;
        go to 120;
      end if
		
			!if the test cannot performed succesfully these are the returned values
			pvalue = 1
			stat = 0
			flag = 0
			
			!remove NA or NULL values from csIndex
			!it can be done in R too

			!if the xIndex is contained in csIndex, x does not bring any new
			!information with respect to cs
			if (ANY(csIndex.eq.xIndex)) then
					pvalues(iter) = 1
                                        stats(iter) = 0;
					go to 120;
			end if

			if(xIndex <= 0) then
					pvalues(iter) = 1
                                        stats(iter) = 0;
					go to 120;
			end if

			!remove dublicates
			!xIndex = unique(xIndex)
			!csIndex = unique(csIndex)

			x = dataset(: , xIndex)

			if(cs_cols>0) then
					cs = dataset(:, csIndex)
			end if

			if(size(x) == 0 .or. size(y) == 0) then
					pvalues(iter) = 1
                                        stats(iter) = 0
					go to 120;
			end if

			!check for constants
			!if(var(x) == 0 .or. var(y) == 0) then
			!        return
			!end if

			!remove constant columns of cs

			if(cs_cols == 0) then!if			
					!calculate the correlation coefficient of x,y directly
					call cor(size(x), x, y, stat);
					stat = abs(stat);
					z = 0.5*log((1+stat)/(1-stat))
					df = size(x) - 3
					w = sqrt(df)*z
					
					call studnt(-abs(w),df, fn_val)
					pvalue = 2*fn_val
			elseif(cs_cols == 1) then
					tmpm(:,1) = x;
					tmpm(:,2) = y;
					tmpm(:,3:cs_cols) = cs;
					call corcoef(tmpm,size(tmpm,1),size(tmpm,2),cor_matrix)

					!perform the test with the cs
					xyIdx = (/ 1,2 /)
					csIdx = (/ 3 /)
					call inv(cor_matrix(csIdx,csIdx), invm, size(invm,1), e)

					rcm = cor_matrix(xyIdx,xyIdx) - MATMUL(cor_matrix(xyIdx,csIdx) , &
					& MATMUL(invm,cor_matrix(csIdx,xyIdx)))

					t1 = rcm(1,2)
					t2 = rcm(1,1)
					t3 = rcm(2,2)

					!stat = abs(rcm(1,2) / sqrt(rcm(1,1) * rcm(2,2)))
					stat = abs(t1 / sqrt(t2 * t3))
					
					!For calculating the stat, I use the below IF because it may be
					!1.0 or a bit above 1 sometimes which is gonna produce a Nan
					!z,w,pvalue due to the log((1+stat)/(1-stat)).
					!SO IN THIS CASE THE RIGHT WAY IS TO USE A COMPLEX TYPE
					!The precision of the digits it may be different from matlab and
					!R but its not a significant difference

					if(stat >= 1.0) then
							!stat = 0.99999
							z_complex = (1+stat)/(1-stat)
							z_complex = 0.5*log(z_complex)
							df = size(x) - cs_cols - 3
							w_complex = sqrt(df)*z_complex

							call studnt(-abs(w_complex),df, fn_val)
							pvalue = 2*fn_val
					else
							z = 0.5*log((1+stat)/(1-stat))
							df = size(x) - cs_cols - 3
							w = sqrt(df)*z
							call studnt(-abs(w),df, fn_val)
							pvalue = 2*fn_val
					endif
			else
					
					tmpm(:,1) = x;
					tmpm(:,2) = y;
					tmpm(:,3:cs_cols) = cs;
					call corcoef(tmpm,size(tmpm,1),size(tmpm,2),cor_matrix)

					xyIdx = (/ 1,2 /)
					!createing the csIdx matrix
					do i=1,cs_cols
							csIdx(i) = i+2;
					end do
					call inv(cor_matrix(csIdx,csIdx), invm, size(invm,1), e)

					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

					rcm = cor_matrix(xyIdx,xyIdx) - MATMUL(cor_matrix(xyIdx,csIdx) ,&
					& MATMUL(invm,cor_matrix(csIdx,xyIdx)))

					t1 = rcm(1,2)
					t2 = rcm(1,1)
					t3 = rcm(2,2)

					!stat = abs(rcm(1,2) / sqrt(rcm(1,1) * rcm(2,2)))
					stat = abs(t1 / sqrt(t2 * t3))

					!For calculating the stat, I use the below IF because it may be
					!1.0 or a bit above 1 sometimes which is gonna produce a Nan
					!z,w,pvalue due to the log((1+stat)/(1-stat)).
					!SO IN THIS CASE THE RIGHT WAY IS TO USE A COMPLEX TYPE
					!The precision of the digits it may be different from matlab and
					!R but its not a significant difference

					if(stat >= 1.0) then
							!stat = 0.99999
							z_complex = (1+stat)/(1-stat)
							z_complex = 0.5*log(z_complex)
							df = size(x) - cs_cols - 3
							w_complex = sqrt(df)*z_complex
							call studnt(-abs(w_complex),df, fn_val)
							pvalue = 2*fn_val
							!pvalue = 2*studnt(-abs(w_complex),df)
					else
							z = 0.5*log((1+stat)/(1-stat))
							df = size(x) - cs_cols - 3
							w = sqrt(df)*z
							call studnt(-abs(w),df, fn_val)
							pvalue = 2*fn_val
					endif

			end if

			flag = 1

			!By definition, NAN is not equal to anything, even itself. Simply compare the variable to itself
			if(pvalue /= pvalue) then
					pvalue = 1;
					stat = 0;
					flag = 0;
			end if
			
			
	120		pvalues(iter) = pvalue
                        stats(iter) = stat
		
        end do

        return

end subroutine fisher_uv
