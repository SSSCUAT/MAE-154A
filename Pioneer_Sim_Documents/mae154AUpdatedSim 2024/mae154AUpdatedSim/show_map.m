%map display function -- D.Toohey
function  out = show_map(in)

pn = in(1);
pe = in(2);
alt = -in(3);
tar_E = in(4);
tar_N = in(5);
whale_E = in(6);
whale_N = in(7);



figure(1)
%subplot(211)

hold on
%hold off
plot(pe,pn,'.b')
%hold on
plot(tar_E,tar_N,'.r','markersize',16);
plot(whale_E,whale_N,'.g','markersize',20);
xlabel('pE')
ylabel('pN')
ylim([-15000 15000]) %used to be [-15000 15000]
xlim([-15000 15000])
axis equal
grid on
legend({'UAV position','Target (whale) position','Whale location'}, 'Location', 'northwest')

% subplot(212)
% hold on
% plot(pn,alt,'.r')
% xlabel('pN (ft)')
% ylabel('Alt (ft)')
% xlim([-15000 15000])
% ylim([500 1500])
% grid on


out = pn;