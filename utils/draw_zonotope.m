function handle = draw_zonotope(G, h, varargin)
Parser = inputParser;
Parser.FunctionName = 'draw_zonotope';
Parser.addOptional('FaceColor', [0.8, 0.2, 0.2]);
Parser.addOptional('FaceAlpha', 0.3);
Parser.parse(varargin{:});

        Vertices = get_all_zonotope_points(G) + h;
        
        indices_convhull = convhull(Vertices');
        Vertices = Vertices(:, indices_convhull);
        
        handle = fill(Vertices(1, :)', Vertices(2, :)', ...
            Parser.Results.FaceColor, ...
            'FaceAlpha', Parser.Results.FaceAlpha); hold on;
end