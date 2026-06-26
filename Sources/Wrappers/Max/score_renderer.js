autowatch = 1;

var scriptVersion = "score_renderer v3";

mgraphics.init();
mgraphics.relative_coords = 0;
mgraphics.autofill = 0;

var svgFilename = "score-max.svg";
var svg = null;
var currentEventNumber = 0;
var maxEventNumber = 255;
var zoomAmount = 1.0;
var offsetX = 0;
var offsetY = 0;
var margin = 10;
var highlightColor = [1, 0, 0, 1];
var blackColor = [0, 0, 0, 1];

post(scriptVersion + " loaded\n");
load_svg();

function load_svg()
{
    svg = new MGraphicsSVG(svgFilename);

    if (!svg || !svg.loaded) {
        post("score_renderer: could not load " + svgFilename + "\n");
    }
}

function paint()
{
    clear_background();

    if (!svg || !svg.loaded) {
        return;
    }

    apply_event_colors();

    var width = box.rect[2] - box.rect[0];
    var height = box.rect[3] - box.rect[1];
    var svgWidth = svg.size[0];
    var svgHeight = svg.size[1];

    if (width <= 0 || height <= 0 || svgWidth <= 0 || svgHeight <= 0) {
        return;
    }

    var availableWidth = Math.max(1, width - (margin * 2));
    var availableHeight = Math.max(1, height - (margin * 2));
    var fit = Math.min(availableWidth / svgWidth, availableHeight / svgHeight);
    var scale = fit * zoomAmount;
    var drawWidth = svgWidth * scale;
    var drawHeight = svgHeight * scale;
    var x = ((width - drawWidth) * 0.5) + offsetX;
    var y = ((height - drawHeight) * 0.5) + offsetY;

    mgraphics.save();
    mgraphics.translate(x, y);
    mgraphics.scale(scale, scale);
    mgraphics.svg_render(svg);
    mgraphics.restore();
}

function apply_event_colors()
{
    svg.mapreset();

    for (var eventNumber = 1; eventNumber <= maxEventNumber; eventNumber++) {
        var targetColor = eventNumber === currentEventNumber ? highlightColor : blackColor;
        svg.mapcolor(event_source_color(eventNumber), targetColor);
    }
}

function event_source_color(eventNumber)
{
    return [0, 0, eventNumber / 255, 1];
}

function clear_background()
{
    mgraphics.set_source_rgba(1, 1, 1, 1);
    mgraphics.paint();
}

function bang()
{
    mgraphics.redraw();
}

function read(filename)
{
    if (filename) {
        svgFilename = filename.toString();
    }

    currentEventNumber = 0;
    load_svg();
    mgraphics.redraw();
}

function msg_int(value)
{
    event(value);
}

function msg_float(value)
{
    event(value);
}

function event(value)
{
    currentEventNumber = Math.max(0, parseInt(value, 10) || 0);
    mgraphics.redraw();
}

function clear()
{
    event(0);
}

function color(red, green, blue)
{
    highlightColor = [
        clamp_color(red) / 255,
        clamp_color(green) / 255,
        clamp_color(blue) / 255,
        1
    ];

    mgraphics.redraw();
}

function clamp_color(value)
{
    return Math.max(0, Math.min(255, parseInt(value, 10) || 0));
}

function zoom(value)
{
    zoomAmount = Math.max(0.01, Number(value) || 1.0);
    mgraphics.redraw();
}

function offset(x, y)
{
    offsetX = Number(x) || 0;
    offsetY = Number(y) || 0;
    mgraphics.redraw();
}

function setmargin(value)
{
    margin = Math.max(0, Number(value) || 0);
    mgraphics.redraw();
}

function onresize()
{
    mgraphics.redraw();
}
