from flask import Flask
app = Flask(__name__)

@app.route('/')
def homepage():
    return """
    <h1>Hello world!</h1>

    <iframe src="http://localhost:8050" width="1253" height="780" frameborder="0" allowfullscreen></iframe>
    """

if __name__ == '__main__':
    app.run(debug=True, use_reloader=True)
