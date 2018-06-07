from remus.remus import app

app.config.update(
    TESTING=True,
    DEGBUG=True,
    ENV="development",
    TRAP_HTTP_EXCEPTIONS=True,
    # USE_X_SENDFILE=True,  # For Future Apache efficient file sending
    EXPLAIN_TEMPLATE_LOADING=True,

)

app.run(debug=True,
        port=5001)
        #threaded=True)
