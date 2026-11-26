// Requires Fiji Bio-Formats. One series, two channels, one Z and T.
path=File.openDialog("Choose two-channel CZI");
run("Bio-Formats Importer", "open=["+path+"] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
getDimensions(w,h,c,z,t);
if(c!=2 || z!=1 || t!=1) exit("Expected C=2, Z=1, T=1");
Dialog.create("Channel mapping");
Dialog.addChoice("P channel",newArray("1","2"),"1"); Dialog.show();
pc=parseInt(Dialog.getChoice()); sc=3-pc;
source=getImageID(); n=w*h;
p=newArray(n); s=newArray(n);
Stack.setChannel(pc); for(i=0;i<n;i++) p[i]=getPixel(i%w,floor(i/w));
Stack.setChannel(sc); for(i=0;i<n;i++) s[i]=getPixel(i%w,floor(i/w));
for(mode=0;mode<2;mode++) {
 selectImage(source); Stack.setChannel(pc); run("Duplicate...", "title=P_test");
 if(mode==0) run("Conversions...", " "); else run("Conversions...", "scale calibrate");
 run("16-bit"); pp=newArray(n); for(i=0;i<n;i++) pp[i]=getPixel(i%w,floor(i/w)); close();
 selectImage(source); Stack.setChannel(sc); run("Duplicate...", "title=S_test"); run("16-bit");
 cp=0; cs=0; cr=0; support=0; ep=0; es=0; er=0; ps=0; ss=0; pn=0; sn=0;
 for(i=0;i<n;i++) {
  q=getPixel(i%w,floor(i/w));
  if(pp[i]!=p[i]) cp++; if(q!=s[i]) cs++;
  ep=maxOf(ep,abs(pp[i]-p[i])); es=maxOf(es,abs(q-s[i]));
  old=p[i]>0 && s[i]>0; now=pp[i]>0 && q>0;
  if(old!=now) support++;
  if(old && now) { e=abs(s[i]/p[i]-q/pp[i]); if(e>1e-10) cr++; er=maxOf(er,e); }
  if(old) {ps+=p[i]; ss+=s[i]; pn+=pp[i]; sn+=q;}
 }
 print("Scaling mode="+mode+" changed P="+cp+" S="+cs+" max errors="+ep+", "+es+" support changes="+support+" ratio changes="+cr+" max ratio error="+er);
 if(ps>0 && pn>0) print("Integrated S/P on original shared support: "+ss/ps+" -> "+sn/pn);
 close();
}
run("Conversions...", " ");
selectImage(source);
