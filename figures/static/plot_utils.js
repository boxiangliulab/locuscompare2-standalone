;(function (root, factory) {
	if (typeof define === 'function' && define.amd) {
		define(['postal'], function () {
			return (root.LocusZoomPlot = factory(root.LocusZoom))
		})
	} else if (typeof module === 'object' && module.exports) {
		module.exports = root.LocusZoomPlot = factory(require('LocusZoom'))
	} else {
		root.LocusZoomPlot = factory(root.LocusZoom)
	}
})(this, function (lz) {
	// Genome base pairs static data
	//    var genome_data = [
	//        { chr: 1, base_pairs: 249250621 },
	//        { chr: 2, base_pairs: 243199373 },
	//        { chr: 3, base_pairs: 198022430 },
	//        { chr: 4, base_pairs: 191154276 },
	//        { chr: 5, base_pairs: 180915260 },
	//        { chr: 6, base_pairs: 171115067 },
	//        { chr: 7, base_pairs: 159138663 },
	//        { chr: 8, base_pairs: 146364022 },
	//        { chr: 9, base_pairs: 141213431 },
	//        { chr: 10, base_pairs: 135534747 },
	//        { chr: 11, base_pairs: 135006516 },
	//        { chr: 12, base_pairs: 133851895 },
	//        { chr: 13, base_pairs: 115169878 },
	//        { chr: 14, base_pairs: 107349540 },
	//        { chr: 15, base_pairs: 102531392 },
	//        { chr: 16, base_pairs: 90354753 },
	//        { chr: 17, base_pairs: 81195210 },
	//        { chr: 18, base_pairs: 78077248 },
	//        { chr: 19, base_pairs: 59128983 },
	//        { chr: 20, base_pairs: 63025520 },
	//        { chr: 21, base_pairs: 48129895 },
	//        { chr: 22, base_pairs: 51304566 },
	//        { chr: 23, base_pairs: 155270560 },
	//        { chr: 24, base_pairs: 59373566 },
	//    ]

	var genome_data = [
		{ chr: 1, base_pairs: 248956422 },
		{ chr: 2, base_pairs: 242193529 },
		{ chr: 3, base_pairs: 198295559 },
		{ chr: 4, base_pairs: 190214555 },
		{ chr: 5, base_pairs: 181538259 },
		{ chr: 6, base_pairs: 170805979 },
		{ chr: 7, base_pairs: 159345973 },
		{ chr: 8, base_pairs: 145138636 },
		{ chr: 9, base_pairs: 138394717 },
		{ chr: 10, base_pairs: 133797422 },
		{ chr: 11, base_pairs: 135086622 },
		{ chr: 12, base_pairs: 133275309 },
		{ chr: 13, base_pairs: 114364328 },
		{ chr: 14, base_pairs: 107043718 },
		{ chr: 15, base_pairs: 101991189 },
		{ chr: 16, base_pairs: 90338345 },
		{ chr: 17, base_pairs: 83257441 },
		{ chr: 18, base_pairs: 80373285 },
		{ chr: 19, base_pairs: 58617616 },
		{ chr: 20, base_pairs: 64444167 },
		{ chr: 21, base_pairs: 46709983 },
		{ chr: 22, base_pairs: 50818468 },
		{ chr: 23, base_pairs: 156040895 },
		{ chr: 24, base_pairs: 57227415 },
	]

	var genome_data_merged = []
	var genome_end = 0
	genome_data.forEach(function (d, i) {
		genome_data_merged.push({})
		genome_data_merged[i].chr = d.chr
		genome_data_merged[i].base_pairs = d.base_pairs

		genome_data_merged[i].genome_start = genome_end

		genome_end += d.base_pairs
		genome_data_merged[i].genome_end = genome_end

		genome_data_merged[i].tickpoint =
			genome_data_merged[i].genome_start + Math.round(d.base_pairs / 2)
	})

	genome_x_list = [
		{
			x: genome_data_merged[0].tickpoint,
			text: '1',
			style: {
				fill: 'rgb(120, 120, 186)',
				'text-anchor': 'center',
				'font-size': '10px',
				'font-weight': 'bold',
			},
			transform: 'translate(0, 5)  rotate(0)',
		},
		{
			x: genome_data_merged[1].tickpoint,
			text: '2',
			style: {
				fill: 'rgb(0, 0, 66)',
				'text-anchor': 'center',
				'font-size': '10px',
				'font-weight': 'bold',
			},
			transform: 'translate(0, 5)  rotate(0)',
		},
		{
			x: genome_data_merged[2].tickpoint,
			text: '3',
			style: {
				fill: 'rgb(120, 120, 186)',
				'text-anchor': 'center',
				'font-size': '10px',
				'font-weight': 'bold',
			},
			transform: 'translate(0, 5)  rotate(0)',
		},
		{
			x: genome_data_merged[3].tickpoint,
			text: '4',
			style: {
				fill: 'rgb(0, 0, 66)',
				'text-anchor': 'center',
				'font-size': '10px',
				'font-weight': 'bold',
			},
			transform: 'translate(0, 5)  rotate(0)',
		},
		{
			x: genome_data_merged[4].tickpoint,
			text: '5',
			style: {
				fill: 'rgb(120, 120, 186)',
				'text-anchor': 'center',
				'font-size': '10px',
				'font-weight': 'bold',
			},
			transform: 'translate(0, 5)  rotate(0)',
		},
		{
			x: genome_data_merged[5].tickpoint,
			text: '6',
			style: {
				fill: 'rgb(0, 0, 66)',
				'text-anchor': 'center',
				'font-size': '10px',
				'font-weight': 'bold',
			},
			transform: 'translate(0, 5)  rotate(0)',
		},
		{
			x: genome_data_merged[6].tickpoint,
			text: '7',
			style: {
				fill: 'rgb(120, 120, 186)',
				'text-anchor': 'center',
				'font-size': '10px',
				'font-weight': 'bold',
			},
			transform: 'translate(0, 5)  rotate(0)',
		},
		{
			x: genome_data_merged[7].tickpoint,
			text: '8',
			style: {
				fill: 'rgb(0, 0, 66)',
				'text-anchor': 'center',
				'font-size': '10px',
				'font-weight': 'bold',
			},
			transform: 'translate(0, 5)  rotate(0)',
		},
		{
			x: genome_data_merged[8].tickpoint,
			text: '9',
			style: {
				fill: 'rgb(120, 120, 186)',
				'text-anchor': 'center',
				'font-size': '10px',
				'font-weight': 'bold',
			},
			transform: 'translate(0, 5)  rotate(0)',
		},
		{
			x: genome_data_merged[9].tickpoint,
			text: '10',
			style: {
				fill: 'rgb(0, 0, 66)',
				'text-anchor': 'center',
				'font-size': '10px',
				'font-weight': 'bold',
			},
			transform: 'translate(0, 5)  rotate(0)',
		},
		{
			x: genome_data_merged[10].tickpoint,
			text: '11',
			style: {
				fill: 'rgb(120, 120, 186)',
				'text-anchor': 'center',
				'font-size': '10px',
				'font-weight': 'bold',
			},
			transform: 'translate(0, 5)  rotate(0)',
		},
		{
			x: genome_data_merged[11].tickpoint,
			text: '12',
			style: {
				fill: 'rgb(0, 0, 66)',
				'text-anchor': 'center',
				'font-size': '10px',
				'font-weight': 'bold',
			},
			transform: 'translate(0, 5)  rotate(0)',
		},
		{
			x: genome_data_merged[12].tickpoint,
			text: '13',
			style: {
				fill: 'rgb(120, 120, 186)',
				'text-anchor': 'center',
				'font-size': '10px',
				'font-weight': 'bold',
			},
			transform: 'translate(0, 5)  rotate(0)',
		},
		{
			x: genome_data_merged[13].tickpoint,
			text: '14',
			style: {
				fill: 'rgb(0, 0, 66)',
				'text-anchor': 'center',
				'font-size': '10px',
				'font-weight': 'bold',
			},
			transform: 'translate(0, 5)  rotate(0)',
		},
		{
			x: genome_data_merged[14].tickpoint,
			text: '15',
			style: {
				fill: 'rgb(120, 120, 186)',
				'text-anchor': 'center',
				'font-size': '10px',
				'font-weight': 'bold',
			},
			transform: 'translate(0, 5)  rotate(0)',
		},
		{
			x: genome_data_merged[15].tickpoint,
			text: '16',
			style: {
				fill: 'rgb(0, 0, 66)',
				'text-anchor': 'center',
				'font-size': '10px',
				'font-weight': 'bold',
			},
			transform: 'translate(0, 5)  rotate(0)',
		},
		{
			x: genome_data_merged[16].tickpoint,
			text: '17',
			style: {
				fill: 'rgb(120, 120, 186)',
				'text-anchor': 'center',
				'font-size': '10px',
				'font-weight': 'bold',
			},
			transform: 'translate(0, 5)  rotate(0)',
		},
		{
			x: genome_data_merged[17].tickpoint,
			text: '18',
			style: {
				fill: 'rgb(0, 0, 66)',
				'text-anchor': 'center',
				'font-size': '10px',
				'font-weight': 'bold',
			},
			transform: 'translate(0, 5)  rotate(0)',
		},
		{
			x: genome_data_merged[18].tickpoint,
			text: '19',
			style: {
				fill: 'rgb(120, 120, 186)',
				'text-anchor': 'center',
				'font-size': '10px',
				'font-weight': 'bold',
			},
			transform: 'translate(0, 5)  rotate(0)',
		},
		{
			x: genome_data_merged[19].tickpoint,
			text: '20',
			style: {
				fill: 'rgb(0, 0, 66)',
				'text-anchor': 'center',
				'font-size': '10px',
				'font-weight': 'bold',
			},
			transform: 'translate(0, 5)  rotate(0)',
		},
		{
			x: genome_data_merged[20].tickpoint,
			text: '21',
			style: {
				fill: 'rgb(120, 120, 186)',
				'text-anchor': 'center',
				'font-size': '10px',
				'font-weight': 'bold',
			},
			transform: 'translate(0, 5)  rotate(0)',
		},
		{
			x: genome_data_merged[21].tickpoint,
			text: '22',
			style: {
				fill: 'rgb(0, 0, 66)',
				'text-anchor': 'center',
				'font-size': '10px',
				'font-weight': 'bold',
			},
			transform: 'translate(0, 5)  rotate(0)',
		},
	]

	lz.Data.EfoWASSource = lz.Data.Source.extend(function (init) {
		this.parseInit(init)
	}, 'EfoWASLZ')

	//we are not using api to fetch data
	lz.Data.EfoWASSource.prototype.getURL = function (state, chain, fields) {
		return ''
	}

	//any post process after getting the JSON
	//we are not using api to fetch data, during init, blank data will be pass to the plot
	//when efo term(s) being selected, these will be replaced with the efo-assoication data
	lz.Data.EfoWASSource.prototype.parseResponse = function (
		resp,
		chain,
		fields,
		outnames,
		trans
	) {
		var data = JSON.parse([])
		return { header: chain.header, body: data }
	}

	//given a location, map it to the plot xais
	transferLocation = function (chr, position) {
		return genome_data_merged[chr - 1].genome_start + position
	}

	reloadlz = function (
		plot_id,
		data_association,
		plotType,
		gene_name,
		colocType,
		profile
	) {
		if (!data_association || data_association.length === 0) {
			console.warn('data is nul')
		}
		let isColoc = plotType === 'coloc'

		//control  x,y,and points
		lz.Layouts.add('data_layer', 'efowas_pvalues', {
			id: 'efowaspvalues',
			type: 'scatter',
			point_shape: {
				scale_function: 'if',
				field: '{{namespace}}category',
				parameters: {
					field_value: 'target_pvalue',
					then: 'diamond',
					else: 'circle',
				},
			},
			point_size: {
				scale_function: 'if',
				field: '{{namespace}}category',
				parameters: {
					field_value: 'target_pvalue',
					then: 80,
					else: isColoc ? 40 : 20,
				},
			},
			tooltip_positioning: 'vertical',
			id_field: '{{namespace}}snp',
			fields: ['{{namespace}}phewas'],
			always_hide_legend: true,
			x_axis: {
				field: '{{namespace}}x',
				floor: 0,
				ceiling: 3094301596,
			},
			y_axis: {
				axis: 1,
				field: '{{namespace}}pval',
				floor: 0,
				upper_buffer: 0.1,
			},
			//xintodo Set defaul color for categories
			color: {
				field: '{{namespace}}category',
				scale_function: 'categorical_bin',
				parameters: {
					null_value: 'rgb(0, 123, 255)',
				},
			},
			legend: [],
			fill_opacity: 1,
			//xintodo datapoint tooltips
			tooltip: {
				closable: true,
				show: { or: ['highlighted', 'selected'] },
				hide: { and: ['unhighlighted', 'unselected'] },
				html:
					'<div> rsid: <strong>{{{{namespace}}snp}}</strong></div>' +
					'<div>pValues: <strong>{{{{namespace}}pvalue}}</strong></div>' +
					'<div>chrom: <strong>{{{{namespace}}chrom}}</strong></div>',
			},
			behaviors: {
				onmouseover: [{ action: 'set', status: 'highlighted' }],
				onmouseout: [{ action: 'unset', status: 'highlighted' }],
				onclick: [
					{ action: 'toggle', status: 'selected', exclusive: true },
				],
				onshiftclick: [{ action: 'toggle', status: 'selected' }],
			},
			label: {
				text: '{{{{namespace}}snp}}',
				spacing: 6,
				lines: {
					style: {
						'stroke-width': '2px',
						stroke: 'red',
						'stroke-dasharray': '20px 20px',
					},
				},
				filters: [
					{
						field: '{{namespace}}pvalue',
						operator: '=',
						value: 7.67957e-9,
					},
				],
				style: {
					'font-size': '14px',
					'font-weight': 'bold',
					fill: '#333333',
				},
			},
		})

		//init x,y ticks, adding data layer
		lz.Layouts.add('panel', 'efowas', {
			id: 'efowas',
			min_width: !isColoc ? 600 : colocType === 'pvalue_plot' ? 500 : 500,
			min_height: !isColoc
				? 200
				: colocType === 'pvalue_plot'
				? 300
				: 150,
			proportional_width: 1,
			// proportional_height: 0.5,
			margin: { top: 0, right: 50, bottom: 40, left: 50 },
			inner_border: 'rgb(210, 210, 210)',
			axes: {
				x: {
					label: isColoc
						? colocType === 'pvalue_plot'
							? profile.trait + ' GWAS -log10(P)'
							: colocType === 'gwas_plot'
							? ''
							: 'Chromosome {{chr}}'
						: 'Genomic Position (number denotes chromosome)',
					label_offset: 35,
					ticks: isColoc ? null : genome_x_list,
					extent: isColoc ? 'state' : null,
				},
				y1: {
					label: !isColoc
						? '-log10 p-value'
						: colocType !== 'gwas_plot'
						? profile.tissue + ' eQTL -log10(P)'
						: profile.trait + ' GWAS -log10(P)',
					label_offset: 28,
				},
			},
			legend: {
				orientation: 'vertical',
				origin: { x: 55, y: 40 },
				hidden: true,
			},
			data_layers: [
				lz.Layouts.get('data_layer', 'significance', {
					unnamespaced: true,
				}),
				lz.Layouts.get('data_layer', 'efowas_pvalues', {
					unnamespaced: true,
				}),
			],
		})

		//the plot, adding efowas panel
		lz.Layouts.add('plot', 'standard_efowas', {
			responsive_resize: false,
			// dashboard: lz.Layouts.get('dashboard', 'standard_plot', {
			//   unnamespaced: true,
			// }),
			panels: [
				lz.Layouts.get('panel', 'efowas', {
					unnamespaced: true,
					proportional_height: 0.45,
				}),
			],
			mouse_guide: false,
		})

		//adding information to association doc
		lz.Data.EfoWASSource.prototype.parseResponse = function (
			resp,
			chain,
			fields,
			outnames,
			trans
		) {
			var data = JSON.parse(JSON.stringify(data_association))
			data = data
				.map((d, i, object) => {
					if (d.chrom) {
						//we only plot those association which have chromosome information
						d.chr = d.chrom
						if (d.chr == 'X') {
							d.chr = 23
						}
						if (d.chr == 'Y') {
							d.chr = 24
						}
						d.bp =
							colocType === 'pvalue_plot'
								? -Math.log10(d.position)
								: parseInt(d.position)
						d.pval = -Math.log10(d.pvalue)
						d.x =
							colocType === 'pvalue_plot'
								? d.bp
								: transferLocation(d.chr, d.bp)
						return d
					} else {
						//some association don't have chromosome information
						// console.warn('no location information available for association:')
						// console.warn(d);
						return false
					}
				})
				.filter((d) => {
					return d
				})
			return { header: chain.header, body: data }
		}

		var data_sources = new lz.DataSources().add('base', [
			'EfoWASLZ',
			{ url: 'www.api.com' },
		])

		var layout = lz.Layouts.get('plot', 'standard_efowas')
		layout.panels[0].margin.top = 10
		layout.panels[0].data_layers[0].offset = 5 // Higher offset for line of GWAS significance than the default 4.522
		if (isColoc) {
			if (colocType === 'pvalue_plot') {
				let chrom = parseInt(data_association[0].chrom)
				layout.state = {
					chr: chrom,
					start: 0,
					end: -Math.log10(data_association[0].position) + 1,
					genome_build: 'GRCh38',
				}
			} else {
				let chrom = parseInt(data_association[0].chrom)
				let start = parseInt(data_association[0].position)
				let end = parseInt(
					data_association[data_association.length - 1].position
				)
				layout.state = {
					chr: chrom,
					start: transferLocation(chrom, start),
					end: transferLocation(chrom, end),
					genome_build: 'GRCh38',
				}
			}

			target_list = data_association.filter(
				(item) => item.category === 'target'
			)
			if (target_list[0] && target_list[0].pvalue) {
				target_pvalue = target_list[0].pvalue
				// 可以显示pvalue>=1e-8的最高点的snp，但是不准，后面有需求再看api
				layout.panels[0].data_layers[1].label.filters[0].value =
					target_pvalue
			}
		}
		layout.panels[0].data_layers[1].color.parameters.categories = [
			'target_pvalue',
			'single_chr',
			'double_chr',
			'r2_five',
			'r2_four',
			'r2_thr',
			'r2_two',
			'r2_one',
		]
		layout.panels[0].data_layers[1].color.parameters.values = [
			'rgb(150, 50, 184)',
			'rgb(175, 175, 175)',
			'rgb(0, 123, 255)',
			'#f50703',
			'#f8a402',
			'#186400',
			'#87ceeb',
			'#0c028b',
		]
		// Generate the plot
		var plot = lz.populate(plot_id, data_sources, layout)
	}

	function manhattanPlot(dom_id, data) {
		console.log('manhattanPlot', data)
		reloadlz(dom_id, data, 'manhattan')
	}

	function colocPlot(dom_id, data, gene_name, plotType, profile) {
		console.log('colocPlot')
		data = data.sort((a, b) => {
			return a.position - b.position
		})
		reloadlz(dom_id, data, 'coloc', gene_name, plotType, profile)
	}

	return { manhattanPlot, colocPlot }
})
