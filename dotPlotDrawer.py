import numpy as np
import plotly.express as px


def drawDotPlot(representationMatrix, location, title):
    img = np.array(representationMatrix) > 0
    fig = px.imshow(img, color_continuous_scale='blues')
    fig.update_xaxes(side='top')
    fig.update_layout(
        title=title,
        width=800,
        height=800,
        coloraxis_showscale=False,
    )
    fig.write_image(location)
