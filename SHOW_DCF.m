for i=1:n
   
    set(h2(i),'XData',pos_x(i)+(range(i)/20)*cos(ph),'YData',pos_y(i)+sin(ph)*(range(i)/20));

    set(h4(i),'XData',pos_x(i)+range(i)*cos(ph),'YData',pos_y(i)+sin(ph)*range(i),'color','c');
   
    if (state(i)==5 | state(i) == 6 | state(i) == 7)

      temp=1-current_frame_length(i)/initial_cfl(i);

      set(h4(i),'XData',pos_x(i)+range(i)*cos(ph),'YData',pos_y(i)+sin(ph)*range(i),'color','m');

      if (within(i,current_frame_dest(i))>0 & state(i)==5)

        set(h5(i),'XData',[pos_x(i),pos_x(i)+min(temp,range(i))*(pos_x(current_frame_dest(i))-pos_x(i))],'YData',[pos_y(i),pos_y(i)+min(temp,range(i))*(pos_y(current_frame_dest(i))-pos_y(i))],'color','r');
        sta_sending = ['Station ', num2str(i) , ' is transmitting to ',num2str(current_frame_dest(i))];
        h1 = legend(sta_sending);
      
      elseif any(state==6)
        set(h5(i),'XData',[pos_x(i),pos_x(i)+min(temp,range(i))*(pos_x(current_frame_dest(i))-pos_x(i))],'YData',[pos_y(i),pos_y(i)+min(temp,range(i))*(pos_y(current_frame_dest(i))-pos_y(i))],'color','m');
        sta_sending = ['Collision occured '];
        h1 = legend(sta_sending);
        
      else
        set(h5(i),'XData',[pos_x(i),pos_x(i)+min(temp,range(i))*(pos_x(current_frame_dest(i))-pos_x(i))],'YData',[pos_y(i),pos_y(i)+min(temp,range(i))*(pos_y(current_frame_dest(i))-pos_y(i))],'color','m');
        sta_sending = ['Out of range '];
        h1 = legend(sta_sending);
      end     
      
   else

      set(h5(i),'color','w');
      
    end 
   
   if ((state(i)==-1 & Timer_SIFS(i)<=0))

       set(h1,'visible','off');
       set(h6(i),'XData',[pos_x(i),pos_x(ACK_dest(i))],'YData',[pos_y(i),pos_y(ACK_dest(i))],'color','g');
        
        sta_sending = ['Station ' num2str(i) ' sending Ack '];
        h1 = legend(sta_sending);
   else
       set(h6(i),'color','w');
   end
   
end %end of for i=1:n


