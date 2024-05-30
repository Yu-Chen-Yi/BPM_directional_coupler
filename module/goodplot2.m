function goodplot2(xlabel_name,ylabel_name,title_name,fontsize)
xlabel(xlabel_name);
ylabel(ylabel_name);
set(gca,'Fontname','times new roman');set(gca,'Fontweight','bold');set(gca,'fontsize',fontsize);
title(title_name,'fontsize',fontsize+4);
end