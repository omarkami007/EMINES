from django.shortcuts import render
import mysql.connector as sql
fn=''
ln=''
em=''
pwd=''
# Create your views here.
def signaction(request):
    global fn,ln,em,pwd
    if request.method=="POST":
        m=sql.connect(host="localhost",user="root",passwd="root",database='test')
        cursor=m.cursor()
        d=request.POST
        for key,value in d.items():
            if key=="prenom":
                fn=value
            if key=="nom":
                ln=value
            if key=="email":
                em=value
            if key=="mot de passe":
                pwd=value
        
        c="insert into users Values('{}','{}','{}','{}')".format(fn,ln,em,pwd)
        cursor.execute(c)
        m.commit()

    return render(request,'signup_page.html')
