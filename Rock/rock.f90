program RockPaperScisors
  implicit none


  integer::counter,reward1,reward2
  real::player1,player2


  counter=0
  do
     counter=counter+1
     call random_number(player1)
     player1=int(3*player1)
     player2=1

     if (player1==0 .and. player2==1) then
        reward1=reward1+1
        reward2=reward2-1
     else if (player1==1 .and. player2==0) then
        reward1=reward1-1
        reward2=reward2+1
     else if (player1==1 .and. player2==2) then
        reward1=reward1+1
        reward2=reward2-1
     else if (player1==2 .and. player2==1) then
        reward1=reward1-1
        reward2=reward2+1
     else if (player1==2 .and. player2==0) then
        reward1=reward1+1
        reward2=reward2-1
     else if (player1==0 .and.player2==2) then
        reward1=reward1-1
        reward2=reward2+1
     end if


     if (mod(counter,10000000)==0)then
        print*, counter
        print*,"Player 1 played",player1
        print*,"Player 1 score:",reward1,reward2,":Player 2 score"
     end if

  end do

end program RockPaperScisors
