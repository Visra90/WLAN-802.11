disp('Length scale of Station spacing is 1');

disp('Please enter the simulation parameters [default value in brackets]...');

n=input('Number of Stations [5]? ');
if length(n)==0
   n=5;
end

r=input('Range of Radio [3]? ');
if length(r)==0
   r=3;
end

motion_scale=input('Time Scale of random motion (SlotTime=2) [100]? ');
if length(motion_scale)==0
   motion_scale=100;
end

frame_size=input('Average Transmission Time (packet size) [8]? ');
if length(frame_size)==0
   frame_size=8;
end

max_simutime=input('Simulation Time in ms (SlotTime=2) [900]? ');%ms
if length(max_simutime)==0
   max_simutime=900;
end

sss=input('Show graphics [y]? ');
if length(sss)==0
   sss='y';
end

disp('Performing calculations -- this may take a while, so please be patient.')
ph=0:0.02:2*pi;

if sss=='y'
figure
hold on
title(' DCF ');
j=sqrt(-1);

    for i=1:n

    h2(i)=plot(exp(j*ph)*(1+r)); %blue node
    h4(i)=plot(exp(j*ph)*(1+r)); %handle for the range ellipse (cyan) - when node is transmitting
    h5(i)=plot(exp(j*ph)*(1+r),'LineWidth',3,'color','r'); %Data packet - RED colour
    h6(i)=plot(exp(j*ph)*(1+r),'LineWidth',3,'color','g'); %ACK packet  - GREEN 
    end
    
end

    alpha=0.5^(1/motion_scale);
    Eb2=0.001;
    beta=(1-alpha^2-Eb2)/(2*Eb2-Eb2*alpha+(1-alpha^2)*alpha);
    ER2=(1-beta^2)*Eb2;

    traffic=1.0;
    ACK_length = 4;
    SIFS=2;
    SlotTime=1;
    DIFS=SIFS+2*SlotTime;
    PCWmin=63;
    PCWmax=1023;
    PCW=(log(PCWmin+1)/log (2))*ones(1,n);

    range=r*ones(1,n);
    state=zeros(1,n);
    prev_state=5*ones(1,n);
    current_frame_length=zeros(1,n);%Sifs+Difs+rand(1,n)*frame_size; 
    Timer_DIFS=zeros(1,n);
    Timer_SIFS=zeros(1,n);
    Timer_ACK=zeros(1,n);
    BC=zeros(1,n);

    pos_x=randn(1,n);
    pos_y=randn(1,n);
    pos_x_change=sqrt(Eb2)*randn(1,n);
    pos_y_change=sqrt(Eb2)*rand(1,n);
    within=ones(n,n);
 
    total_transmissions = 0;
    successful_transmission = 0;
    total_collisions = 0;
    unreachable_packets = 0;
    total_acks = 0;
    ack_collisions = 0;
    unreachable_acks = 0;
    successful_acks = 0;
 
for counter=1:max_simutime
    t0=clock;
    
    pos_x_change=beta*pos_x_change+randn(1,n)*sqrt(ER2);
    pos_y_change=beta*pos_y_change+randn(1,n)*sqrt(ER2);
    pos_x=pos_x*alpha+pos_x_change;
    pos_y=pos_y*alpha+pos_y_change;
    within=ones(n,n);
    
    for i = 1:n
	      for k = 1:n        
		             if((pos_x(i)-pos_x(k))^2 + (pos_y(i)-pos_y(k))^2 > range(i)^2)                     
					                within(i,k) = 0;            
                     end
          end
    end
   
  %% State - 0
   
    temp=rand(1,n);
    transit0to1=zeros(1,n);
    transit0to1(state==0 & temp<traffic)=1;
    
% select back off counter from 0 to minimum contention window

    temp=floor(rand(1,n).*(2.^PCW-1));
    BC(transit0to1>0 & state==0) = temp(transit0to1>0 & state==0)*SlotTime;

%     for j=1:n
%         if BC(j)<=-1
%             BC(j)=temp(j)*SlotTime;
%         end
%     end


%% State - 1
    % Medium Sensing 
    
   sending=zeros(1,n);
   sending(state>=4 | state<=-1 )=1;
   % Busy_media determines which nodes are within range and sending
   Busy_media=sending*within;
   
   transit1to2 = zeros(1,n);     
   transit1to2(state == 1 & Busy_media < 1) = 1;
   
   
%% State - 2     
    % Difs Timer
   
    Timer_DIFS(transit1to2>0) = DIFS+1; % Wait for "DIFS" amount of time before sending the Data packet. DIFS == Distributed coordination inter-frame space
    
    transit2to1 = zeros(1,n);
    transit2to1(state==2 & Busy_media>0)=1;    
    
     % Counting down for DIFS amount of time
    Timer_DIFS(state==2) = Timer_DIFS(state==2)-1;
    
       if (state==2)
       if (Timer_DIFS == DIFS)
        disp('wait DIFS before initial transmit');
       end
        disp(Timer_DIFS(state==2));
       end

     transit2to3 = zeros(1,n);
     transit2to3(state== 2 & Timer_DIFS<0) = 1;
     
 %% State - 3   
   
 % BC -> Backoff Counter
   BC(state==3 & Busy_media<1) = BC(state==3 & Busy_media<1)-1;
   %print the backoff counter
   if (state==3 & Busy_media<1)
      if (Timer_DIFS==0)
        disp('ALL stations Decrement their backoff counter');
      end
   disp(BC(state==3 & Busy_media<1));
   end
   
   transit3to1(state==3 & Busy_media>0)=1;
       
   transit3to5 = zeros(1,n);
   transit3to5(state==3 & BC<0) = 1;
   
   if any(transit3to5>0)
    
        temp=SIFS+rand(1,n)*frame_size*2;
        current_frame_length=temp;
        initial_cfl=current_frame_length;

        temp=ceil(rand(1,n)*n);
        current_frame_dest = temp;

        for i=1:n
           while (current_frame_dest(i)==i | current_frame_dest(i)==0)
              current_frame_dest(i)=ceil(rand(1,1)*n);
           end
        end
   end
  %% State - 5    
   %transmission
  
    for i=1:n   
        for j=1:n
       
            if (state(i)==5 & state(j)==5 & j~=i)
                state(j)=6;
                state(i)=6;
                if prev_state(i)==5
                disp([' OOPS COLLISION OCCURS: station ' num2str(j) ' is also transmitting']);
                else
                disp([' OOPS COLLISION OCCURS: station ' num2str(i) ' is also transmitting']);
                end
            end
        end

       if (state(i)==5 & within(i,current_frame_dest(i))<1)
          state(i)=7;
          disp(['station ' num2str(current_frame_dest(i)) ' is out of ranges from station ' num2str(i)]);
      end

    end

    current_frame_length(state==5 |state==6 | state==7) = current_frame_length(state==5 | state==6 | state==7)-1;

    for i=1:n
        if(state(i) == 5 & current_frame_length(i) > 0)
            disp(['station ' num2str(i) ' sending data to station ' num2str(current_frame_dest(i))])                
        end
    end

     transit5to0 = zeros(1,n);
     transit5to0(state >=5 & current_frame_length < 0) = 1;
     
     transit5to4 = zeros(1,n);
     transit6to0 = zeros(1,n);
     transit7to0 = zeros(1,n);
     
     for i = 1:n
         if(transit5to0(i) > 0)
             if(state(i)==6)
             PCW(state == 6) = PCW(state == 6)+1;
             transit6to0((i)) = 1;
             elseif(state(i)==7)
             transit7to0((i)) = 1;
             else
             transit5to4(current_frame_dest(i)) = 1;
             Timer_SIFS(current_frame_dest(i)) = SIFS+1;  
             Timer_ACK(current_frame_dest(i)) = ACK_length;
             ACK_dest(current_frame_dest(i)) = i;
             end
          end
     end
     
    for cont=1:n
         if (state(cont)==-2 | state(cont)==6)
            cont_size=2.^PCW(cont)-1;
            disp(['contention window of station ' num2str(cont) ' increases to ' num2str(cont_size) ])
        end
    end

    maxpcw=(log(PCWmax+1)/log (2))*ones(1,n);

    if PCW(state == -2 | state ==6)>= maxpcw(state == -2 | state ==6)
        PCW(state == -2 | state ==6)= maxpcw(state == -2 | state ==6);
        disp ('maximum contention window reached');
    end

     total_transmissions = total_transmissions + length(state(state >= 5 & transit5to0 > 0));
     successful_transmission = successful_transmission + length(state(state == 5 & transit5to0 >0));
     total_collisions = total_collisions + length(state(state == 6 & transit5to0 > 0));
     unreachable_packets = unreachable_packets + length(state(state == 7 & transit5to0 > 0));

    state(transit5to4 > 0) = 4;
    Timer_SIFS(state == 4 & transit5to4 > 0) = SIFS+1;     
    Timer_SIFS(state == 4) = Timer_SIFS(state == 4) - 1;
    
    % print the SIFS count down   
    for j=1:n
        if (state(j)==4 & Timer_SIFS(j)>=-1)
            if (Timer_SIFS(j)==SIFS)
            disp(['station ' num2str(j) ' wait SIFS before sending ACK']);
            end
       disp(Timer_SIFS(state==4 & Timer_SIFS(j)>=0));
        end
    end
     
     transit4tom1 = zeros(1,n);
     transit4tom1(state == 4 & Timer_SIFS <= 0) = 1; 
     
  %% State - -1
     
     for i = 1:n
              
         if (state(i) == -1 && Busy_media(ACK_dest(i)) > 1)
             state(i) = -2;
         end
         
         if (state(i) == -1 && within(i,ACK_dest(i)) < 1)                         
             state(i) = -3;
         end
     end
     
    Timer_ACK(state <= -1) = Timer_ACK(state <= -1) - 1;
   
    for j=1:n
        if (state(j)==-1 & Timer_SIFS(j)<=0)
            if Timer_SIFS(j)==0
            disp(['station ' num2str(j) ' sending ACK to ' num2str(ACK_dest(j)) ]);
            end
            disp(Timer_ACK(state==-1 & Timer_SIFS<=0));
        end
    end

    
     transitm1to0 = zeros(1,n);
     transitm1to0(state <= -1 & Timer_ACK < 0) = 1;
     
     total_acks = total_acks + length(state(state <= -1 & transitm1to0 > 0));
     successful_acks = successful_acks + length(state(state == -1 & transitm1to0 > 0));
     ack_collisions = ack_collisions + length(state(state == -2 & transitm1to0 > 0));
     unreachable_acks = unreachable_acks + length(state(state == -3 & transitm1to0 > 0));
  
 %% States Conversion
     
     state(transit0to1 > 0)  =  1;
     state(transit1to2 > 0)  =  2;
     state(transit2to3 > 0)  =  3;
     state(transit3to5 > 0)  =  5;     
     state(transit2to1 > 0)  =  1;
     state(transit5to0 > 0)  =  0;
     state(transit6to0 > 0)  =  0;
     state(transit7to0 > 0)  =  0;
     state(transit5to4 > 0)  =  4;
     state(transit4tom1 > 0) = -1;
     state(transitm1to0 > 0) =  0;
     
     prev_state=state;
     if sss=='y'
     show_DCF
     drawnow

while etime(clock,t0)<0.1
end
     end   
end

 total_transmissions
 successful_transmission
 total_collisions
 unreachable_packets
 total_acks
 successful_acks
 ack_collisions
 unreachable_acks 
   %STATE TABLE 
%-6: Waiting for ACK
%-3: DATA collision?
%-2: ACK collision?
%-1: sending an ACK
%0: idle (nothing to send)
%1: Waiting for free media (have something to send)
%2: Free media detected, wait for DIFS
%3: Random back-off
%4: Waiting for SIFS to send ACK
%5: sending successfully
%6: sending to a busy station
%7: sending to an out-of-range node
%8: sending an RTS signal/packet
%9: sending a CTS signal/packet
   