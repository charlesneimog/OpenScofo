autowatch = 1;

mgraphics.init();
mgraphics.relative_coords = 0;
mgraphics.autofill = 0;

var svgFilename = "score-max.svg";
var svg = null;
var svgWidth = 1;
var svgHeight = 1;
var currentEventNumber = 0;
var maxEventNumber = 255;
var zoomAmount = 1.0;
var offsetX = 0;
var offsetY = 0;
var margin = 10;
var highlightColor = [1, 0, 0, 1];
var blackColor = [0, 0, 0, 1];

load_svg();

function load_svg()
{
    var dimensions = read_svg_dimensions(svgFilename);
    if (!dimensions) {
        svg = null;
        post("score_renderer: could not load " + svgFilename + "\n");
        return;
    }

    svgWidth = dimensions[0];
    svgHeight = dimensions[1];

    try {
        svg = new MGraphicsSVG(svgFilename);
    } catch (error) {
        svg = null;
        post("score_renderer: could not load " + svgFilename + "\n");
    }
}

function read_svg_dimensions(filename)
{
    var file = new File(filename, "read");
    if (!file.isopen) {
        return null;
    }

    var header = file.readstring(Math.min(file.eof, 8192));
    file.close();

    var number = "([+-]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?)";
    var viewBox = new RegExp(
        "viewBox\\s*=\\s*[\\\"']\\s*" + number + "[\\s,]+" + number + "[\\s,]+" + number + "[\\s,]+" + number,
        "i"
    ).exec(header);

    if (viewBox) {
        var viewBoxWidth = Number(viewBox[3]);
        var viewBoxHeight = Number(viewBox[4]);
        if (isFinite(viewBoxWidth) && isFinite(viewBoxHeight) && viewBoxWidth > 0 && viewBoxHeight > 0) {
            return [viewBoxWidth, viewBoxHeight];
        }
    }

    var width = /<svg[\s\S]*?\bwidth\s*=\s*["']\s*([0-9]*\.?[0-9]+)/i.exec(header);
    var height = /<svg[\s\S]*?\bheight\s*=\s*["']\s*([0-9]*\.?[0-9]+)/i.exec(header);
    if (width && height) {
        var documentWidth = Number(width[1]);
        var documentHeight = Number(height[1]);
        if (documentWidth > 0 && documentHeight > 0) {
            return [documentWidth, documentHeight];
        }
    }

    return null;
}

function paint()
{
    clear_background();

    if (!svg) {
        return;
    }

    apply_event_colors();

    var width = box.rect[2] - box.rect[0];
    var height = box.rect[3] - box.rect[1];
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

    mgraphics.svg_render(svg, x, y, drawWidth, drawHeight);
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
