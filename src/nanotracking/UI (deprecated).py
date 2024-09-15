import bottle
from bottle import route, post, redirect, run, template, BaseTemplate, static_file, request
from .DifferencePlotter import *

app = bottle.default_app()
BaseTemplate.defaults['get_url'] = app.get_url

@route('/static/<filepath:path>', name = 'static')
def server_static(filepath):
    return static_file(filepath, root='/Volumes/LAB DRIVE/240829NTA with folder structure/')

@route('/')
def index():
    query_number = request.query.number
    number = 10 if query_number=='' else query_number
    # <p name="number">The number is {{number_val}}</p>
    # <form action="/buttonpress?number={{number_val}}" method="POST">
    #     <td><input type ="submit" name="hello" value="Hello!"></td>
    # </form>
    return template(r'''
    <form action="/rerun_program" method="POST">
        <td><input type ="submit" name="rerun_program" value="Run program"></td>
    </form>
    <img src="/static/{{ img_src }}" style="max-width:400%;max-height:400%;">
    ''', number_val = number, img_src = 'CSV outputs/Ridgeline plot.png')
    # return template(r'''<form action="/buttonpress" method="POST">
    #     <td><input type ="submit" name="Hello!" value="Hello!"></td>
    # </form>''')

@post('/buttonpress')
def buttonpress():
    # hello = request.forms.hello
    # request.number += 1
    # return "<p>Yay!</p>"
    number = int(request.query.number)
    number += 10
    # return f"Number is {number}"
    redirect(f'/?number={number}')

@post('/rerun_program')
def rerun_program():
    run_program()
    # redirect(request.path)
    redirect('/')


if __name__ == '__main__':
    run(host='localhost', port=8080, debug=True, reloader=True)