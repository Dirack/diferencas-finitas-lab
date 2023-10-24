# Makefile
# 
# Objetivo: Compilar e rodar o experimento numérico deste diretório.
# 
# Site: https://dirack.github.io
# 
# Versão 1.0
# 
# Programador: Rodolfo A C Neves (Dirack) 20/10/2023
# 
# Email: rodolfo_profissional@hotmail.com
# 
# Licença: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.


all:	fdmodeling.x experiment view

fdmodeling.x:	Mfdmodeling.c
	gcc -g $< -I../include -L../lib -o $@ -lrsf -lm

experiment:
	../bin/sfmath output=2.0 n1=200 n2=200 d1=0.005 d2=0.005 label1=x1 unit1=km label2=x2 unit2=km --out=stdout > vel.rsf
	< vel.rsf ./fdmodeling.x border=y nb=30 nt=1000 dt=0.001 rec=receiver_with_border.rsf --out=stdout > wav_com_borda.rsf
	< wav_com_borda.rsf ../bin/sfgrey gainpanel=all title="Wavefield Propagation" scalebar=y color=j scalebar=y > wav_com_borda.vpl

view:
	< wav_com_borda.vpl ../bin/sfpen

clean:
	rm *.x *.rsf *.vpl

help:
	@echo "Experimento de modelagem com diferenças finitas"
	@echo "Para compilar o programa de modelagem utilize 'make fdmodeling.x'"
	@echo "Para rodar a modelagem utilize 'make experiment'"
	@echo "Para vizualizar o experimento utilize 'make view'"
	@echo "Para rodar os 3 comandos acima em sequência utilize 'make'"
